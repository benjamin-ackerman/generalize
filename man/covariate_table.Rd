% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariate_table.R
\name{covariate_table}
\alias{covariate_table}
\title{Create Covariate Balance Table}
\usage{
covariate_table(
  trial,
  selection_covariates,
  data,
  weighted_table = FALSE,
  selection_method = "lr",
  sl_library = NULL,
  survey_weights = FALSE,
  trim_weights = FALSE,
  trim_pctile = 0.97,
  is_data_disjoint = TRUE
)
}
\arguments{
\item{trial}{variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)}

\item{selection_covariates}{vector of covariate names in data set that predict trial participation}

\item{data}{data frame comprised of "stacked" trial and target population data}

\item{weighted_table}{should the output be a weighted table?
If \code{TRUE}, then \code{\link{weighting}}}

\item{selection_method}{method to estimate the probability of trial
participation.  Default is logistic regression ("lr").
Other methods supported are Random Forests ("rf") and Lasso ("lasso"),
passed to \code{\link{weighting}}}

\item{sl_library}{vector of SuperLearner library methods. If `selection_method` = 'super', specify names of methods to include in library. Default is NULL.}

\item{survey_weights}{variable name of population data's complex survey weights. Default is \code{FALSE}: if \code{FALSE}, then population data do not come a complex survey and weights do not need to be incorporated in estimation.}

\item{trim_weights}{logical. If \code{TRUE}, then trim the weights to the value specified in `trim_pctile`. Default is \code{FALSE}.}

\item{trim_pctile}{numeric. If `trim_weights` is \code{TRUE}, then specify what percentile weights should be trimmed to. Default is 0.97.}

\item{is_data_disjoint}{logical. If \code{TRUE}, then trial and population data
are considered independent.  This affects calculation of the weights -
see details for more information.}
}
\description{
Create Covariate Balance Table
}
