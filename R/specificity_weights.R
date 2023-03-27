


#' Calculates specificity scores for a gene (from a genes x cell-type matrix) using SPM method
#' as described in Xiao SJ, Zhang C, Zou Q, Ji ZL. TiSGeD: a database for tissue-specific genes. Bioinformatics. 2010 May 1;26(9):1273-5. doi: 10.1093/bioinformatics/btq109. Epub 2010 Mar 11. PMID: 20223836; PMCID: PMC2859128.
#' code adapted from Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform. 2017 Mar 1;18(2):205-214. doi: 10.1093/bib/bbw008. PMID: 26891983; PMCID: PMC5444245.
#'
#' @param data a vector of gene expression values across cell-types
#'
#' @return a specificity score value (double)
#'
#' @export
#'
Spm <- function(data){
  if(all(!is.na(data)))
  {
    if(min(data, na.rm=TRUE) >= 0)
    {
      if(sum(data) !=0)
      {
        spm <- data^2/(data%*%data)
        res <- max(spm)
      } else {
        res <- 0
      }
    } else {
      res <- NA
    }
  } else {
    res <- NA
  }
  return(res)
}


#' Calculates specificity scores for a gene (from a genes x cell-type matrix) using Hg method (entropy)
#' as described in Schug J, Schuller WP, Kappen C, Salbaum JM, Bucan M, Stoeckert CJ Jr. Promoter features related to tissue specificity as measured by Shannon entropy. Genome Biol. 2005;6(4):R33. doi: 10.1186/gb-2005-6-4-r33. Epub 2005 Mar 29. PMID: 15833120; PMCID: PMC1088961.
#' code adapted from Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform. 2017 Mar 1;18(2):205-214. doi: 10.1093/bib/bbw008. PMID: 26891983; PMCID: PMC5444245.
#'
#' @param data a vector of gene expression values across cell-types
#'
#' @return a specificity score value (double)
#'
#' @export
#'
Hg <- function(data){
  if(all(!is.na(data)))
  {
    if(min(data, na.rm=TRUE) >= 0)
    {
      if(sum(data) !=0)
      {
        p <- data / sum(data)
        res <- -sum(p*log2(p), na.rm=TRUE)
        res <- 1 - (res/log2(length(p))) #Modification: To bring to normalized scale
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}


#' calculate specificity scores
#'
#' @param dataset a genes x cell-types matrix from the reference data
#' @param method the specificity weight calculation method used, currently SPM and Hg methods are supported
#' @param transform an option for transforming the data into log or sqrt space prior to calculating specificity scores
#'
#' @return a vector of specificity weights for each gene
#'
#' @export
#'

calc_specificity_weights <- function(dataset, method = 'SPM', transform = 'log'){
  if(!(transform %in% c('none', 'log', 'sqrt'))){
    stop('current transformation options are none, log, and sqrt')
  }
  if(!(method %in% c('SPM', 'Hg'))){
    stop('current specificity weighting method options are SPM and Hg')
  }

  if(transform == 'log'){
    if(max(dataset) < 50){warning('has the dataset already been log transformed? If so, set transform = "none"')}
    dataset <- log1p(dataset)
  }
  if(transform == 'sqrt'){
    if(max(dataset) < 250){warning('has the dataset already been sqrt transformed? If so, set transform = "none"')}
    dataset <- sqrt(dataset)
  }

  if(method == 'SPM'){
    weights <- pbapply::pbapply(dataset, 1, Spm)
  }
  if(method == 'Hg'){
    weights <- pbapply::pbapply(dataset, 1, Hg)
  }

  return(weights)
}

