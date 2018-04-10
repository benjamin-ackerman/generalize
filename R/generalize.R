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
#' @examples


generalize <- function(outcome, treatment, trial, selection_covariates, data, method = "weighting",
                       selection_method = "lr", is_data_disjoint = TRUE, trim_pop = TRUE){

  ##### make methods lower case #####
  method = tolower(method)
  selection_method = tolower(selection_method)

  ##### CHECKS #####
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  if(anyNA(match(outcome,names(data)))){
    stop("Outcome is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(treatment,names(data)))){
    stop("Treatment is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(unique(data[,trial])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(anyNA(match(names(table(data[,trial])),c("0","1")))){
    stop("Sample Membership variable must be coded as `0` (not in trial) or `1` (in trial)",call. = FALSE)
  }

  if(!length(unique(data[,treatment])) == 2){
    stop("Treatment variable not binary", call. = FALSE)
  }

  if(!method %in% c("weighting","bart","tmle")){
    stop("Invalid method!",call. = FALSE)
  }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid weighting method!",call. = FALSE)
  }

  ##### just keep the data we need #####
  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]

  ##### trim population #####
  if(trim_pop == TRUE){
    data = trim_pop(trial, covariates, data)
  }

  ##### estimate SATE #####
  SATE_model = lm(as.formula(paste(outcome,treatment,sep="~")), data = data)

  SATE = summary(SATE_model)$coefficients[treatment, "Estimate"]
  SATE_se = summary(SATE_model)$coefficients[treatment, "Std. Error"]

  SATE_CI_l = SATE - 1.96*SATE_se
  SATE_CI_u = SATE + 1.96*SATE_se

  SATE_results = c(SATE,SATE_se,SATE_CI_l,SATE_CI_u)

  ##### Weighting Methods #####
  if(method == "weighting"){
    TATE_results = weighting(outcome, treatment, trial, selection_covariates, data, selection_method, is_data_disjoint, trim_pop)$TATE
  }

  if(method == "BART"){
    TATE_results = "NOT READY YET"
  }

  if(method == "TMLE"){
    TATE_results = tmle(outcome, treatment, trial, selection_covariates, data)$TATE
  }

  result.tab = rbind(SATE_results, TATE_results)
  colnames(result.tab) = c("Estimate","Std. Error","95% CI Lower","95% CI Upper")
  row.names(result.tab) = c("SATE","TATE")

  return(result.tab)
}
