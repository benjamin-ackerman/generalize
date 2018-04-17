
#' Generalize Average Treatment Effect from Randomized Trial to Population
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param method method to generalize average treatment effect to the target population.  Default is "weighting" (weighting by participation probability).  Other methods supported are "BART" (Bayesian Additive Regression Trees - NOT READY YET) and "TMLE" (Targeted Maximum Likelihood Estimation)
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param trim_pop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of trial covariates.
#' @return \code{generalize} returns an object of the class "generalize"
# @examples
# generalize(outcome = "studyretention", treatment = "treat", trial = "trial", selection_covariates = c("age","sex","race","ethnicity","marital_status"),data = meth_data)
# generalize(outcome = "studyretention", treatment = "treat", trial = "trial", selection_covariates = c("age","sex","race","ethnicity","marital_status"),data = meth_data,method = "tmle")
# summary(generalize(outcome = "studyretention", treatment = "treat", trial = "trial", selection_covariates = c("age","sex","race","ethnicity","marital_status"), data = meth_data, method = "weighting", selection_method = "rf", trim_pop = TRUE))

generalize <- function(outcome, treatment, trial, selection_covariates, data, method = "weighting",
                       selection_method = "lr", is_data_disjoint = TRUE, trim_pop = FALSE){

  ##### make methods lower case #####
  method = tolower(method)
  selection_method = tolower(selection_method)

  ### If using TMLE, must trim target population
  if(method == "tmle"){
    trim_pop = TRUE
  }

  ##### CHECKS #####
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  if(is.null(outcome) | anyNA(match(outcome,names(data)))){
    stop("Outcome is not a variable in the data provided!",call. = FALSE)
  }

  if(is.null(treatment) | anyNA(match(treatment,names(data)))){
    stop("Treatment is not a variable in the data provided!",call. = FALSE)
  }

  if(is.null(selection_covariates) | anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(na.omit(unique(data[,trial]))) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(anyNA(match(names(table(data[,trial])),c("0","1")))){
    stop("Sample Membership variable must be coded as `0` (not in trial) or `1` (in trial)",call. = FALSE)
  }

  if(!length(na.omit(unique(data[,treatment]))) == 2){
    stop("Treatment variable not binary", call. = FALSE)
  }

  if(!method %in% c("weighting","bart","tmle")){
    stop("Invalid method!",call. = FALSE)
  }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid weighting method!",call. = FALSE)
  }

  ##### trim population #####
  if(trim_pop == FALSE){
    n_excluded = NULL
    ## just keep the data we need
    data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]
  }

  if(trim_pop == TRUE){
    n_excluded = trim_pop(trial, selection_covariates, data)$n_excluded
    data = trim_pop(trial, selection_covariates, data)$trimmed_data
  }

  ##### Weighting object for diagnostics #####
  weight_object = weighting(outcome, treatment, trial, selection_covariates, data, selection_method, is_data_disjoint)

  participation_probs = weight_object$participation_probs
  weights = weight_object$weights
  g_index = gen_index(participation_probs$population, participation_probs$trial)

  ##### Generalize results #####
  ## First, estimate SATE
  SATE_model = lm(as.formula(paste(outcome,treatment,sep="~")), data = data)

  SATE = summary(SATE_model)$coefficients[treatment, "Estimate"]
  SATE_se = summary(SATE_model)$coefficients[treatment, "Std. Error"]

  SATE_CI_l = SATE - 1.96*SATE_se
  SATE_CI_u = SATE + 1.96*SATE_se

  SATE_results = list(estimate = SATE,
                      se = SATE_se,
                      CI_l = SATE_CI_l,
                      CI_u = SATE_CI_u)

  ## Weighting results
  if(method == "weighting"){
    TATE_results = weight_object$TATE
  }

  ## BART results
  if(method == "bart"){
    TATE_results = "NOT READY YET"
  }

  ## TMLE results
  if(method == "tmle"){
    TATE_results = tmle(outcome, treatment, trial, selection_covariates, data)$TATE
  }

  ##### sample size of trial and population #####
  n_trial = nrow(data[which(data[,trial] == 1),])
  n_pop = nrow(data[which(data[,trial] == 0),])

  #cov_tab = covariate_table(trial, selection_covariates, data)

  data_output = data[,c(outcome, treatment, trial, selection_covariates)]

  ##### Items to save to "generalize" object #####
  out = list(
    SATE = SATE_results,
    TATE = TATE_results,
    outcome = outcome,
    treatment = treatment,
    trial = trial,
    method = method,
    selection_method = selection_method,
    g_index = g_index,
    n_trial = n_trial,
    n_pop = n_pop,
    trim_pop = trim_pop,
    n_excluded = n_excluded,
    selection_covariates = selection_covariates,
    #covariate_table = cov_tab,
    data = data_output,
    is_data_disjoint = is_data_disjoint
  )

  class(out) = "generalize"

  return(out)
}

print.generalize <- function(x,...){
  cat("A generalize object: \n")
  cat(paste0(" - SATE: ", round(x$SATE$estimate,3), "\n"))
  cat(paste0(" - TATE: ", round(x$TATE$estimate,3), "\n"))
  cat(paste0(" - outcome variable: ", x$outcome, "\n"))
  cat(paste0(" - treatment variable: ", x$treatment, "\n"))
  cat(paste0(" - generalizability method: ", x$method, "\n"))
  if(x$method == "weighting"){
    cat(paste0("     - probability of trial participation method: ", x$selection_method, "\n"))
  }
  cat(paste0(" - common covariates included: ", paste(x$selection_covariates, collapse = ", "), "\n"))
  cat(paste0(" - sample size of trial: ", x$n_trial, "\n"))
  cat(paste0(" - size of population: ", x$n_pop, "\n"))
  cat(paste0(" - was population trimmed according to trial covariate bounds?: ", ifelse(x$trim_pop == TRUE, "Yes", "No"), "\n"))
  if(x$trim_pop == TRUE){
    cat(paste0("    - number excluded from population data: ", x$n_excluded, "\n"))
  }

  invisible(x)
}

summary.generalize <- function(object,...){
  ## put together results table
  result_tab = rbind(unlist(object$SATE), unlist(object$TATE))
  colnames(result_tab) = c("Estimate","Std. Error","95% CI Lower","95% CI Upper")
  row.names(result_tab) = c("SATE","TATE")

  ## give full names to methods
  method_name = c("Weighting", "TMLE", "BART")
  method = c("weighting","tmle","bart")

  selection_method_name = c("Logistic Regression","Random Forests","Lasso")
  selection_method = c("lr","rf","lasso")

  ## build outcome formula
  outcome_formula = paste0(object$outcome, " ~ ", object$treatment)

  out = list(
    outcome_formula = outcome_formula,
    result_tab = result_tab,
    method = method_name[object$method == method],
    selection_method = selection_method_name[object$selection_method == selection_method],
    n_trial = object$n_trial,
    n_pop = object$n_pop,
    trim_pop = object$trim_pop,
    n_excluded = object$n_excluded,
    g_index = object$g_index
  )

  class(out) = "summary.generalize"

  return(out)
}

print.summary.generalize <- function(x,...){
  cat("Average Treatment Effect Estimates: \n \n")
  cat(paste0("Outcome Model: ",x$outcome_formula," \n \n"))
  print(x$result_tab)

  cat("\n")
  cat("============================================ \n")
  cat(paste0("TATE estimated by ",x$method, "\n"))
  if(x$method == "Weighting"){
    cat(paste0("Weights estimated by ", x$selection_method,"\n"))
  }
  cat("\n")
  cat(paste0("Trial sample size: ",x$n_trial,"\n"))
  cat(paste0("Population size: ",x$n_pop,"\n"))
  if(x$trim_pop == TRUE){
    cat("Population data were trimmed for covariates to not exceed trial covariate bounds \n")
    cat(paste0("Number excluded from population: ", x$n_excluded ,"\n"))
  }
  cat("\n")
  cat(paste0("Generalizability Index: ", round(x$g_index,3), "\n"))

  #cat("Covariate Distributions: \n \n")
  #print(round(x$covariate_table,4))

  invisible(x)
}

