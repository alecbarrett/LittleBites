% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subtraction.R
\name{get_best_LR_single_log}
\alias{get_best_LR_single_log}
\title{calculates the ROC-AUC value for a single bulk sample}
\usage{
get_best_LR_single_log(
  bulk_vector,
  sample_name,
  pseudobulk,
  cells_in_use,
  training_matrix,
  training_genes,
  anti_identity_vector,
  proportions_table,
  LR_list,
  specificity_score,
  sample_name_separator = "r",
  fast_auc = TRUE
)
}
\arguments{
\item{bulk_vector}{a vector of gene expression in the sample (with named entries)}

\item{sample_name}{a character vector, the name of a sample being processed}

\item{pseudobulk}{a genes x cell-types matrix of gene expression values from the single cell reference}

\item{cells_in_use}{a character vector of cell types being used for proportion modeling and subtraction}

\item{training_matrix}{a genes x cell-types ground truth matrix}

\item{training_genes}{a vector of gene names to calculate the AUC with}

\item{anti_identity_vector}{a vector of ones and zeros that matches the proportions vector in naming and length, with a zero value that corresponds to the target cell type (so that the target cell type is not subtracted)}

\item{proportions_table}{a table of estimated proportions for the bulk sample (with entries named for the cell-types)}

\item{LR_list}{list of values to use as the learning rates: ex: 1/(2**seq(0,10,1)) returns a list of 1, 1/2, 1/4, 1/8, etc...}

\item{specificity_score}{a vector of gene level weights calculated using the specificity score from the single cell reference}

\item{sample_name_separator}{a character value to split the cell type from the replicate number in the sample name (assumed structure is cell-sample_name_separator-number)}

\item{fast_auc}{if true, then function will use a new faster way to calculate auc than the original implementation}
}
\value{
a vector of gene expression values, named by the genes
}
\description{
calculates the ROC-AUC value for a single bulk sample
}
