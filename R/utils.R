
utils::globalVariables(c("threshold"))

#' Calculates TPR when given a prediction matrix and a ground truth matrix
#'
#' @param expression genes x cell-types gene expression matrix
#' @param truth genes x cell-types ground truth matrix
#' @param threshold a value to set the prediction for the expression matrix
#'
#' @return a TPR value for the full matrix
#'
#' @export
#'
get_tpr <- function(expression,
                    truth,
                    threshold){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  if(length(expression) != length(truth)){
    stop('expression and truth vectors must be the same length')
  }
  bin <- expression >= threshold

  return(sum(bin * truth)/sum(truth))
}


#' Calculates FPR when given a prediction matrix and a ground truth matrix
#'
#' @param expression genes x cell-types gene expression matrix
#' @param truth genes x cell-types ground truth matrix
#' @param threshold a value to set the prediction for the expression matrix
#'
#' @return an FPR value for the full matrix
#'
#' @export
#'
get_fpr <- function(expression,
                    truth,
                    threshold){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  if(length(expression) != length(truth)){
    stop('expression and truth vectors must be the same length')
  }

  bin <- expression >= threshold

  return(sum(bin * (!truth))/sum(!(truth)))
}


#' Calculates FDR when given a prediction matrix and a ground truth matrix
#'
#' @param expression genes x cell-types gene expression matrix
#' @param truth genes x cell-types ground truth matrix
#' @param threshold a value to set the prediction for the expression matrix
#'
#' @return an FDR value for the full matrix
#'
#' @export
#'
get_fdr <- function(expression,
                    truth,
                    threshold){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  if(length(expression) != length(truth)){
    stop('expression and truth vectors must be the same length')
  }

  bin <- expression >= threshold
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))

  return(fdr)
}


#' calculates the ROC-AUC value for a single bulk sample
#'
#' @param bulk_vector a vector of gene expression in the sample (with named entries)
#' @param sample_name a character vector, the name of a sample being processed
#' @param training_matrix a genes x cell-types ground truth matrix
#' @param training_genes a vector of gene names to calculate the AUC with
#' @param threshold_list a vector of values to threshold the gene expression dataset, used to calculate the ROC-AUC
#' @param sep a character value to split the cell type from the replicate number in the sample name (assumed structure is cell-sep-number)
#'
#' @return ROC-AUC value for this bulk sample
#'
#' @export
#'
calc_bulk_auc_one_sample <- function(bulk_vector,
                                     sample_name,
                                     training_matrix,
                                     training_genes,
                                     threshold_list = c(0,seq(0,15,0.1)),
                                     sep = 'r'){

  # log transform the bulk data, then calculate the AUC for the ROC curve.

  ## bulk_vector = genes x samples matrix, needs to be named with gene names
  ## training_matrix = genes x bulk_cell_types matrix
  ## training_genes = list of genes

  cell_type <- stringr::str_split_fixed(sample_name, sep, 2)[,1]
  dcnt <- log1p(bulk_vector)
  names(dcnt) <- names(bulk_vector)
  dcnt <- dcnt[training_genes]


  train_bulk <- training_matrix[training_genes,]


  train_bulk_single <- train_bulk[,cell_type]
  names(train_bulk_single) <- rownames(train_bulk)
  diags_d <- tibble::tibble(threshold = threshold_list,
                    TPR = purrr::map_dbl(threshold, ~get_tpr(dcnt, train_bulk_single, .x)),
                    FPR = purrr::map_dbl(threshold, ~get_fpr(dcnt, train_bulk_single, .x)))

  bulk_auc <- bayestestR::auc(rev(diags_d$FPR), rev(diags_d$TPR) , method = 'trap')
  return(bulk_auc)
}


# fast auc calculation functions for littlebites algorithm
# these functions provide significant speedup while maintaining identical results

#' vectorized tpr calculation for multiple thresholds
#'
#' @param expression vector of expression values
#' @param truth vector of ground truth (0/1 or logical)
#' @param thresholds vector of threshold values
#'
#' @return vector of tpr values for each threshold
#'
#' @export
#'
get_tpr_vectorized <- function(expression, truth, thresholds) {
  # convert to matrix: genes x thresholds
  n_genes <- length(expression)
  n_thresholds <- length(thresholds)
  
  # create matrices for vectorized comparison
  expr_matrix <- matrix(expression, nrow = n_genes, ncol = n_thresholds)
  thresh_matrix <- matrix(thresholds, nrow = n_genes, ncol = n_thresholds, byrow = TRUE)
  
  # vectorized binary classification
  bin_matrix <- expr_matrix >= thresh_matrix
  
  # vectorized tpr calculation
  P <- sum(truth)
  if (P == 0) return(rep(0, n_thresholds))
  
  TPR <- colSums(bin_matrix * truth) / P
  return(TPR)
}

#' vectorized fpr calculation for multiple thresholds
#'
#' @param expression vector of expression values
#' @param truth vector of ground truth (0/1 or logical)
#' @param thresholds vector of threshold values
#'
#' @return vector of fpr values for each threshold
#'
#' @export
#'
get_fpr_vectorized <- function(expression, truth, thresholds) {
  # convert to matrix: genes x thresholds
  n_genes <- length(expression)
  n_thresholds <- length(thresholds)
  
  # create matrices for vectorized comparison
  expr_matrix <- matrix(expression, nrow = n_genes, ncol = n_thresholds)
  thresh_matrix <- matrix(thresholds, nrow = n_genes, ncol = n_thresholds, byrow = TRUE)
  
  # vectorized binary classification
  bin_matrix <- expr_matrix >= thresh_matrix
  
  # vectorized fpr calculation
  N <- sum(!truth)
  if (N == 0) return(rep(0, n_thresholds))
  
  FPR <- colSums(bin_matrix * (!truth)) / N
  return(FPR)
}

#' fast auc calculation using vectorized threshold processing
#'
#' @param expression vector of expression values
#' @param truth vector of ground truth (0/1 or logical)
#' @param threshold_list vector of thresholds to use
#'
#' @return auc value
#'
#' @export
#'
calc_auc_fast <- function(expression, truth, threshold_list = c(0, seq(0, 15, 0.1))) {
  # handle edge cases
  if (length(expression) != length(truth)) {
    stop('expression and truth vectors must be the same length')
  }
  
  if (sum(truth) == 0 || sum(truth) == length(truth)) {
    return(0.5)  # no discrimination possible
  }
  
  # vectorized tpr/fpr calculation
  TPR <- get_tpr_vectorized(expression, truth, threshold_list)
  FPR <- get_fpr_vectorized(expression, truth, threshold_list)
  
  # use bayestestR::auc exactly like the original
  auc <- bayestestR::auc(rev(FPR), rev(TPR), method = 'trap')
  
  return(auc)
}

#' fast bulk auc calculation for a single sample
#'
#' @param bulk_vector a vector of gene expression in the sample (with named entries)
#' @param sample_name a character vector, the name of a sample being processed
#' @param training_matrix a genes x cell-types ground truth matrix
#' @param training_genes a vector of gene names to calculate the auc with
#' @param threshold_list a vector of values to threshold the gene expression dataset
#' @param sep a character value to split the cell type from the replicate number
#'
#' @return roc-auc value for this bulk sample
#'
#' @export
#'
calc_bulk_auc_one_sample_fast <- function(bulk_vector,
                                         sample_name,
                                         training_matrix,
                                         training_genes,
                                         threshold_list = c(0, seq(0, 15, 0.1)),
                                         sep = 'r') {
  
  # extract cell type from sample name
  cell_type <- stringr::str_split_fixed(sample_name, sep, 2)[,1]
  
  # log transform and subset
  dcnt <- log1p(bulk_vector)
  names(dcnt) <- names(bulk_vector)
  dcnt <- dcnt[training_genes]
  
  # get ground truth
  train_bulk <- training_matrix[training_genes,]
  train_bulk_single <- train_bulk[,cell_type]
  names(train_bulk_single) <- rownames(train_bulk)
  
  # use fast vectorized auc calculation
  bulk_auc <- calc_auc_fast(dcnt, train_bulk_single, threshold_list)
  
  return(bulk_auc)
}




