% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specificity_weights.R
\name{calc_specificity_weights}
\alias{calc_specificity_weights}
\title{calculate specificity scores}
\usage{
calc_specificity_weights(dataset, method = "SPM", transform = "log")
}
\arguments{
\item{dataset}{a genes x cell-types matrix from the reference data}

\item{method}{the specificity weight calculation method used, currently SPM and Hg methods are supported}

\item{transform}{an option for transforming the data into log or sqrt space prior to calculating specificity scores}
}
\value{
a vector of specificity weights for each gene
}
\description{
calculate specificity scores
}
