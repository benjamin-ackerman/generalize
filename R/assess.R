#' Assess Generalizability of Randomized Trial to Population
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param trim_pop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of trial covariates.
#' @param seed numeric. By default, the seed is set to 13783, otherwise can be specified (such as for simulation purposes).
#' @return \code{assess} returns an object of the class "generalize_assess"

assess = function(trial, selection_covariates, data, selection_method = "lr",
                  is_data_disjoint = TRUE, trim_pop = FALSE,seed){

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

  if(!missing(seed)){
    if(!is.numeric(seed)){
      stop("seed must be numeric!,call. = FALSE")
    }}
  ##### Clean up data from missing values #####
  #data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),]

  if(trim_pop == FALSE){n_excluded = NULL}

  if(trim_pop == TRUE){
    n_excluded = trim_pop(trial, selection_covariates, data)$n_excluded
    data = trim_pop(trial, selection_covariates, data)$trimmed_data
  }

  weight_object = weighting(outcome = NULL, treatment = NULL, trial, selection_covariates, data,
                               selection_method, is_data_disjoint,seed)

  participation_probs = weight_object$participation_probs
  weights = weight_object$weights

  g_index = gen_index(participation_probs$trial, participation_probs$population)

  cov_tab = covariate_table(trial, selection_covariates, data)
  weighted_cov_tab = covariate_table(trial, selection_covariates, data,weighted_table = TRUE)


  n_trial = nrow(data[which(data[,trial] == 1),])
  n_pop = nrow(data[which(data[,trial] == 0),])

  data_output = data[,c(trial, selection_covariates)]

  out = list(
    g_index = g_index,
    selection_method = selection_method,
    selection_covariates = selection_covariates,
    trial_name = trial,
    n_trial = n_trial,
    n_pop = n_pop,
    trim_pop = trim_pop,
    n_excluded = n_excluded,
    participation_probs = participation_probs,
    weights = weights,
    covariate_table = cov_tab,
    weighted_covariate_table = weighted_cov_tab,
    data = data_output
    )

  class(out) = "generalize_assess"

  return(out)
}

print.generalize_assess <- function(x,...){
  cat("A generalize_assess object: \n")
  cat(paste0(" - probability of trial participation method: ", x$selection_method, "\n"))
  cat(paste0(" - common covariates included: ", paste(x$selection_covariates, collapse = ", "), "\n"))
  cat(paste0(" - sample size of trial: ", x$n_trial, "\n"))
  cat(paste0(" - size of population: ", x$n_pop, "\n"))
  cat(paste0(" - was population trimmed according to trial covariate bounds?: ", ifelse(x$trim_pop == TRUE, "Yes", "No"), "\n"))
  if(x$trim_pop == TRUE){
    cat(paste0("    - number excluded from population data: ", x$n_excluded, "\n"))
    }

  invisible(x)
}

summary.generalize_assess <- function(object,...){
  selection_method_name = c("Logistic Regression","Random Forests","Lasso")
  selection_method = c("lr","rf","lasso")
  prob_dist_table = rbind(summary(object$participation_probs$trial),
                          summary(object$participation_probs$population))
  row.names(prob_dist_table) = paste0(c("Trial","Population"), " (n = ", c(object$n_trial,object$n_pop),")")

  selection_formula = paste0(object$trial_name," ~ ",paste(object$selection_covariates, collapse = " + "))

  out = list(
    selection_formula = selection_formula,
    selection_method = selection_method_name[selection_method == object$selection_method],
    g_index = round(object$g_index,3),
    prob_dist_table = prob_dist_table,
    covariate_table = round(object$covariate_table, 4),
    weighted_covariate_table = round(object$weighted_covariate_table,4),
    trim_pop = object$trim_pop,
    n_excluded = object$n_excluded
  )

  class(out) = "summary.generalize_assess"
  return(out)
}

print.summary.generalize_assess <- function(x,...){
  cat("Probability of Trial Participation: \n \n")
  cat(paste0("Selection Model: ",x$selection_formula," \n \n"))
  print(x$prob_dist_table)
  cat("\n")
  cat(paste0("Estimated by ",x$selection_method, "\n"))
  cat(paste0("Generalizability Index: ", round(x$g_index,3), "\n"))
  cat("============================================ \n")
  cat("Covariate Distributions: \n \n")
  if(x$trim_pop == TRUE){
    cat("Population data were trimmed for covariates to not exceed trial covariate bounds \n")
    cat(paste0("Number excluded from population: ", x$n_excluded ,"\n \n"))
  }
  print(round(x$covariate_table,4))
  invisible(x)
}
