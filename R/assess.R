#' Assess Generalizability of Randomized Trial to Population
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param sl_library vector of SuperLearner library methods. If `selection_method` = 'super', specify names of methods to include in library. Default is NULL.
#' @param survey_weights variable name of population data's complex survey weights. Default is FALSE: if FALSE, then population data do not come a complex survey and weights do not need to be incorporated in estimation.
#' @param trim_weights logical. If TRUE, then trim the weights to the value specified in `trim_pctile`. Default is FALSE.
#' @param trim_pctile numeric. If `trim_weights` is TRUE, then specify what percentile weights should be trimmed to. Default is 0.97.
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param trimpop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of trial covariates.
#' @param seed numeric. By default, the seed is set to 13783, otherwise can be specified (such as for simulation purposes).
#' @return \code{assess} returns an object of the class "generalize_assess"

#' @export
assess = function(trial, selection_covariates, data, selection_method = "lr",
                  sl_library = NULL, survey_weights = FALSE, trim_weights=FALSE, trim_pctile = .97,
                  is_data_disjoint = TRUE, trimpop = FALSE,seed){

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

  if(!selection_method %in% c("lr","rf","lasso","gbm","super")){
    stop("Invalid weighting method!",call. = FALSE)
  }

  if(!missing(seed)){
    if(!is.numeric(seed)){
      stop("seed must be numeric!,call. = FALSE")
    }}
  ##### Clean up data from missing values #####
  #data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),]

  if(trimpop == FALSE){n_excluded = NULL}

  if(trimpop == TRUE){
    n_excluded = trim_pop(trial, selection_covariates, data)$n_excluded
    data = trim_pop(trial, selection_covariates, data)$trimmed_data
  }

  weight_object = weighting(outcome=NULL, treatment=NULL, trial, selection_covariates, data, selection_method, sl_library,
                            survey_weights, trim_weights, trim_pctile, is_data_disjoint,seed)

  participation_probs = weight_object$participation_probs
  weights = weight_object$weights

  g_index = gen_index(participation_probs$trial, participation_probs$population)

  cov_tab = covariate_table(trial = trial, selection_covariates = selection_covariates, data = data,
                                     weighted_table = FALSE, survey_weights=survey_weights)

  weighted_cov_tab = covariate_table(trial = trial, selection_covariates = selection_covariates, data = data,
                                     weighted_table = TRUE, selection_method = selection_method, sl_library = sl_library, survey_weights=survey_weights,
                                     trim_weights=trim_weights, trim_pctile=trim_pctile,is_data_disjoint = is_data_disjoint)


  n_trial = nrow(data[which(data[,trial] == 1),])
  n_pop = nrow(data[which(data[,trial] == 0),])

  if(survey_weights == FALSE){
    data_output = data[,c(trial, selection_covariates)]
    n_pop_eff = NA
  } else{
    data_output = data[,c(trial, selection_covariates,survey_weights)]
    n_pop_eff = ceiling(sum(data[which(data[,trial]==0),survey_weights], na.rm=TRUE))
  }


  out = list(
    g_index = g_index,
    selection_method = selection_method,
    selection_covariates = selection_covariates,
    trial_name = trial,
    n_trial = n_trial,
    n_pop = n_pop,
    n_pop_eff = n_pop_eff,
    trimpop = trimpop,
    n_excluded = n_excluded,
    participation_probs = participation_probs,
    weights = weights,
    survey_weights = (survey_weights != FALSE),
    covariate_table = cov_tab,
    weighted_covariate_table = weighted_cov_tab,
    data = data_output
    )

  class(out) = "generalize_assess"

  return(out)
}

#' @export
print.generalize_assess <- function(x,...){
  cat("A generalize_assess object: \n")
  cat(paste0(" - probability of trial participation method: ", x$selection_method, "\n"))
  cat(paste0(" - common covariates included: ", paste(x$selection_covariates, collapse = ", "), "\n"))
  cat(paste0(" - sample size of trial: ", x$n_trial, "\n"))
  cat(paste0(" - size of population: ", x$n_pop, "\n"))
  if(x$survey_weights == TRUE){
    cat(paste0(" - size of effective population (accounting for survey weights): ",x$n_pop_eff,"\n"))
  }
  cat(paste0(" - was population trimmed according to trial covariate bounds?: ", ifelse(x$trimpop == TRUE, "Yes", "No"), "\n"))
  if(x$trimpop == TRUE){
    cat(paste0("    - number excluded from population data: ", x$n_excluded, "\n"))
    }

  invisible(x)
}

#' @export
summary.generalize_assess <- function(object,...){
  selection_method_name = c("Logistic Regression","Random Forests","Lasso","GBM","SuperLearner")
  selection_method = c("lr","rf","lasso","gbm","super")

  prob_dist_table = rbind(summary(object$participation_probs$trial),
                          summary(object$participation_probs$population))
  row.names(prob_dist_table) = paste0(c("Trial","Population"), " (n = ", c(object$n_trial,object$n_pop),")")

  selection_formula = paste0(object$trial_name," ~ ",paste(object$selection_covariates, collapse = " + "))

  out = list(
    selection_formula = selection_formula,
    selection_method = selection_method_name[selection_method == object$selection_method],
    survey_weights = object$survey_weights,
    n_pop_eff = object$n_pop_eff,
    g_index = round(object$g_index,3),
    prob_dist_table = prob_dist_table,
    covariate_table = round(object$covariate_table, 4),
    weighted_covariate_table = round(object$weighted_covariate_table,4),
    trimpop = object$trimpop,
    n_excluded = object$n_excluded
  )

  class(out) = "summary.generalize_assess"
  return(out)
}

#' @export
print.summary.generalize_assess <- function(x,...){
  cat("Probability of Trial Participation: \n \n")
  cat(paste0("Selection Model: ",x$selection_formula," \n \n"))
  print(x$prob_dist_table)
  cat("\n")
  if(x$survey_weights == TRUE){
    cat(paste0(" - Effective population size (accounting for survey weights): ",x$n_pop_eff,"\n"))
    cat("\n")
  }
  cat(paste0("Estimated by ",x$selection_method, "\n"))
  cat(paste0("Generalizability Index: ", round(x$g_index,3), "\n"))
  cat("============================================ \n")
  cat("Covariate Distributions: \n \n")
  if(x$trimpop == TRUE){
    cat("Population data were trimmed for covariates to not exceed trial covariate bounds \n")
    cat(paste0("Number excluded from population: ", x$n_excluded ,"\n \n"))
  }
  print(round(x$covariate_table,4))
  invisible(x)
}
