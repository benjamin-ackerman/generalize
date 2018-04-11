#' Assess Generalizability of Randomized Trial to Population
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param trim_pop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of trial covariates.
#' @return \code{generalize} returns an object of the class "generalize"
#' @examples

assess = function(trial, selection_covariates, data, selection_method = "lr",
                  is_data_disjoint = TRUE, trim_pop = FALSE){

  ##### make methods lower case #####
  selection_method = tolower(selection_method)

  ##### CHECKS #####
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  if(anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(na.omit(unique(data[,trial]))) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(anyNA(match(names(table(data[,trial])),c("0","1")))){
    stop("Sample Membership variable must be coded as `0` (not in trial) or `1` (in trial)",call. = FALSE)
  }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid weighting method!",call. = FALSE)
  }

  ##### Clean up data from missing values #####
  #data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),]

  if(trim_pop == FALSE){n_excluded = NULL}

  if(trim_pop == TRUE){
    n_excluded = trim_pop(trial, selection_covariates, data)$n_excluded
  }

  weighting_object = weighting(outcome = NULL, treatment = NULL, trial, selection_covariates, data,
                               selection_method, is_data_disjoint, trim_pop)

  participation_probs = weighting_object$participation_probs
  weights = weighting_object$weights

  g_index = gen_index(participation_probs$probs_trial, participation_probs$probs_population)

  out = list(
    n_excluded = n_excluded,
    g_index = g_index)

  return(out)
}
