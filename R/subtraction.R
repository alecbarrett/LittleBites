#' Function to iteratively subtract across samples
#'
#' @param bulk a vector of gene expression in the sample (with named entries)
#' @param reference a genes x cell-types matrix of gene expression values from the single cell reference
#' @param sep a character value to split the cell type from the replicate number in the sample name (assumed structure is cell-sep-number)
#' @param cell_types_matrix a samples x cell-types matrix describing the target and contaminant cell types to estimate per sample
#' @param training_matrix a genes x cell-types ground truth matrix
#' @param specificity_weights a vector of gene level weights calculated using the specificity score from the single cell reference
#' @param learning_rate list of values to use as the learning rates: ex: 1/(2**seq(0,10,1)) returns a list of 1, 1/2, 1/4, 1/8, etc...
#' @param max_iterations maximum number of iterations to go through
#' @param verbose if true, then the function will return information about each iteration
#'
#' @return a matrix of cleaned bulk expression
#'
#' @export
#'


subtraction <- function(bulk,
                        reference,
                        sep = 'r',
                        cell_types_matrix,
                        training_matrix,
                        specificity_weights,
                        learning_rate = c(0, seq(1/(2**seq(0,10,1)))),
                        max_iterations = 100,
                        verbose = T){

  if(sep == ''){stop('function requires a separator')}

  bulk_deconv <- bulk

  samples <- colnames(bulk)

  for(sample_1 in samples){
    if(verbose){print(sample_1)}

    cells_in_use <- cell_types_matrix[sample_1,] |> unlist()

    cell <- cells_in_use[1]
    contaminant_tissues <- cells_in_use[2:length(cells_in_use)]

    for(i in seq(0,max_iterations,1)){

      if(i==0){

        bulk_deconv_target <- bulk[,sample_1] ## some steps required a dataframe
        names(bulk_deconv_target) <- rownames(bulk)
        starting_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                                 sample_1,
                                                 training_matrix,
                                                 training_genes = rownames(training_matrix),

                                                 sep = sep)
        pre_auc <- starting_auc

      if(verbose){print(c('starting AUC', starting_auc))}
      }


      else{
        if(verbose){print(c('iteration',i))}

        reference_tmp <- reference[,cells_in_use]
        reference_tmp <- sweep(reference_tmp, 2, colSums(reference_tmp), '/')
        reference_tmp <- reference_tmp * sum(bulk_deconv_target)
        estimates <- nnls::nnls( A = log1p(as.matrix(reference_tmp * specificity_weights )),
                           b = log1p(bulk_deconv_target * specificity_weights ) )$x
        names(estimates) <- cells_in_use
        estimates <- estimates/sum(estimates)
        if(verbose){print(estimates)}


        specificity_weights_use <- specificity_weights[rownames(reference_tmp)]


        anti_identity_vector <- (names(estimates) %in% contaminant_tissues) * 1

        bulk_deconv_target <- get_best_LR_single_log(bulk_deconv_target,
                                                     sample_1,
                                                     reference_tmp,
                                                     cells_in_use,
                                                     training_matrix,
                                                     training_genes = rownames(training_matrix),
                                                     anti_identity_vector,
                                                     data.frame(estimates),
                                                     40/(2**seq(0,20,1)),
                                                     specificity_weights_use,
                                                     sep = sep)


        bulk_deconv_target <- bulk_deconv_target[[1]]
        if(verbose){print(c('subtracted ', calc_bulk_auc_one_sample(bulk_deconv_target,
                                                                    sample_1,
                                                                    training_matrix,
                                                                    training_genes = rownames(training_matrix),
                                                                    sep = sep)))}




        post_auc <- calc_bulk_auc_one_sample(bulk_deconv_target,
                                             sample_1,
                                             training_matrix,
                                             training_genes = rownames(training_matrix),
                                             sep = sep)

        if(verbose){print(c(i, post_auc))}


        if(post_auc == pre_auc){
          if(verbose){print(c('total AUC percent improvement -->',(post_auc - starting_auc)*100))}
          bulk_deconv[,sample_1] <- bulk_deconv_target

          reference_tmp <- reference[,cells_in_use]
          reference_tmp <- sweep(reference_tmp, 2, colSums(reference_tmp), '/')
          reference_tmp <- reference_tmp * sum(bulk_deconv_target)
          estimates <- nnls::nnls( A=log1p(as.matrix(reference_tmp * specificity_weights )),
                             b=log1p(bulk_deconv_target * specificity_weights ) )$x
          names(estimates) <- cells_in_use
          estimates <- estimates/sum(estimates)
          if(verbose){print('final proportion estimates:')}
          if(verbose){print(estimates)}

          break}
        pre_auc <- post_auc
      }
    }
  }
  print('done')
  return(bulk_deconv)
}



#' log transforms bulk and reference profiles, and the subtracts contaminating profiles
#'
#' @param learning_rate a value used to scale the subtraction magnitude
#' @param bulk a vector of gene expression in the sample (with named entries)
#' @param pseudobulk a genes x cell-types matrix of gene expression values from the single cell reference
#' @param cells_in_use a character vector of cell types being used for proportion modeling and subtraction
#' @param proportions_vector a vector of estimated proportions for the bulk sample (with entries named for the cell-types)
#' @param anti_identity_vector a vector of ones and zeros that matches the proportions vector in naming and length, with a zero value that corresponds to the target cell type (so that the target cell type is not subtracted)
#' @param specificity_score a vector of gene level weights calculated using the specificity score from the single cell reference
#'
#' @return a vector of gene expression values, named by the genes
#'
#' @export
#'
subtract_single_log <- function(learning_rate,
                                bulk,
                                pseudobulk,
                                cells_in_use,
                                proportions_vector,
                                anti_identity_vector,
                                specificity_score){

  ## learning_rate = scalar from 0 to 1
  ## bulk = genes x samples matrix
  ## pseudbulk = genes x sc_cell_types matrix
  ## proportions_table = sc_cell_types x samples matrix
  #pseudobulk <- sc_pseudobulk_use
  sc.log <- log(pseudobulk[,cells_in_use])
  sc.log[sc.log<0]=0
  sc.log <- sc.log * specificity_score

  proportions_vector <- proportions_vector * anti_identity_vector

  subtraction_vector <- (learning_rate * (as.matrix(sc.log) %*%  as.matrix(proportions_vector)))

  bulk_deconv <- log(bulk) - subtraction_vector
  bulk_deconv <- exp(bulk_deconv)

  return(bulk_deconv)
}

#' calculates the ROC-AUC value for a single bulk sample
#'
#' @param bulk_vector a vector of gene expression in the sample (with named entries)
#' @param sample_name a character vector, the name of a sample being processed
#' @param pseudobulk a genes x cell-types matrix of gene expression values from the single cell reference
#' @param cells_in_use a character vector of cell types being used for proportion modeling and subtraction
#' @param training_matrix a genes x cell-types ground truth matrix
#' @param training_genes a vector of gene names to calculate the AUC with
#' @param anti_identity_vector a vector of ones and zeros that matches the proportions vector in naming and length, with a zero value that corresponds to the target cell type (so that the target cell type is not subtracted)
#' @param proportions_table a table of estimated proportions for the bulk sample (with entries named for the cell-types)
#' @param LR_list list of values to use as the learning rates: ex: 1/(2**seq(0,10,1)) returns a list of 1, 1/2, 1/4, 1/8, etc...
#' @param specificity_score a vector of gene level weights calculated using the specificity score from the single cell reference
#' @param sep a character value to split the cell type from the replicate number in the sample name (assumed structure is cell-sep-number)
#'
#' @return a vector of gene expression values, named by the genes
#'
#' @export
#'
get_best_LR_single_log <- function(bulk_vector,
                                   sample_name,
                                   pseudobulk,
                                   cells_in_use,
                                   training_matrix,
                                   training_genes,
                                   anti_identity_vector,
                                   proportions_table,
                                   LR_list,
                                   specificity_score,
                                   sep = 'r'){

  # take a list of learning rates, and find the best one as defined by the one with the highest AUC value
  ## bulk = vector of gene values for a single sample
  ## pseudbulk = genes x sc_cell_types matrix
  ## training_matrix = genes x bulk_cell_types matrix
  ## training_genes = list of genes
  ## proportions table = sc_cell_types estimates for a single sample
  ## LR_list = list of values to use as the learning rates: ex: 1/(2**seq(0,10,1)) returns a list of 1, 1/2, 1/4, 1/8, etc...




  if((0 %in% LR_list) == F){LR_list <- c(0,LR_list)}
  bulk_deconv_list <- lapply(LR_list, function(LR){
    #print(LR)
    bulk_subtract_df <- subtract_single_log(learning_rate = LR,
                                            bulk = bulk_vector,
                                            pseudobulk = pseudobulk,
                                            cells_in_use = cells_in_use,
                                            proportions_vector = proportions_table,
                                            anti_identity_vector = anti_identity_vector,
                                            specificity_score = specificity_score)
    bulk_subtract <- bulk_subtract_df[,1]
    names(bulk_subtract) <- rownames(bulk_subtract_df)
    bulk_subtract <- bulk_subtract[order(names(bulk_subtract))]
    return(bulk_subtract)

  })

  auc_list <- lapply(bulk_deconv_list, function(bulk_subtracted){ ## for all elements in bulk_deconv_list, calculate their AUC
    auc_add <- calc_bulk_auc_one_sample(bulk_subtracted,
                                        sample_name,
                                        training_matrix,
                                        training_genes = rownames(training_matrix),
                                        sep = sep)
    #print(auc_add)
    return(auc_add)
  })
  best_LR <- which(unlist(auc_list) == max(unlist(auc_list)), arr.ind = TRUE)
  #print(best_LR)
  if(1 %in% best_LR){
    best_LR <- 1
  }
  else { best_LR <- utils::tail(best_LR, n = 1) }
  return(list(bulk_deconv_list[best_LR][[1]], LR_list[best_LR])) ## return a matrix with the highest AUC
  #return(bulk_deconv_list)
}
