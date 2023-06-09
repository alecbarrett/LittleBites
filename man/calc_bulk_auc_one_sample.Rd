% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_bulk_auc_one_sample}
\alias{calc_bulk_auc_one_sample}
\title{calculates the ROC-AUC value for a single bulk sample}
\usage{
calc_bulk_auc_one_sample(
  bulk_vector,
  sample_name,
  training_matrix,
  training_genes,
  threshold_list = c(0, seq(0, 15, 0.1)),
  sep = "r"
)
}
\arguments{
\item{bulk_vector}{a vector of gene expression in the sample (with named entries)}

\item{sample_name}{a character vector, the name of a sample being processed}

\item{training_matrix}{a genes x cell-types ground truth matrix}

\item{training_genes}{a vector of gene names to calculate the AUC with}

\item{threshold_list}{a vector of values to threshold the gene expression dataset, used to calculate the ROC-AUC}

\item{sep}{a character value to split the cell type from the replicate number in the sample name (assumed structure is cell-sep-number)}
}
\value{
ROC-AUC value for this bulk sample
}
\description{
calculates the ROC-AUC value for a single bulk sample
}
