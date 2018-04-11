<<<<<<< HEAD
##### MAYBE USE THIS FOR EVERYTHING #####
## Can output table of TATE results
## Can generate participation probs and weights to be used for
    ## generalizability index
    ## plotting the distributions
    ## KEEP ASSESS
### ***write code to make outcome and treatment null, then just use for assessing generalizability***

=======
>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
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
<<<<<<< HEAD
=======
#' generalize(outcome = "STUDYCOMPLETE", treatment = "treat", trial = "trial", selection_covariates = c("age","sex","race"), data = ctn_data, method = "weighting", selection_method = "rf")

>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4

generalize <- function(outcome, treatment, trial, selection_covariates, data, method = "weighting",
                       selection_method = "lr", is_data_disjoint = TRUE, trim_pop = TRUE){

<<<<<<< HEAD
generalize <- function(outcome, treatment, trial, selection_covariates, data, method = "weighting",
                       selection_method = "lr", is_data_disjoint = TRUE, trim_pop = FALSE){

=======
>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
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

<<<<<<< HEAD

  if(anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(na.omit(unique(data[,trial]))) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(anyNA(match(names(table(data[,trial])),c("0","1")))){
    stop("Sample Membership variable must be coded as `0` (not in trial) or `1` (in trial)",call. = FALSE)
  }

  if(length(na.omit(unique(data[,treatment]))) == 2){
=======
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
>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
    stop("Treatment variable not binary", call. = FALSE)
  }

  if(!method %in% c("weighting","bart","tmle")){
    stop("Invalid method!",call. = FALSE)
  }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid weighting method!",call. = FALSE)
  }
<<<<<<< HEAD


  ##### trim population #####
  if(trim_pop == FALSE){n_excluded = NULL}

  if(trim_pop == TRUE){
    n_excluded = trim_pop(trial, selection_covariates, data)$n_excluded
    data = trim_pop(trial, selection_covariates, data)$trimmed_data
  }

  ##### just keep the data we need #####
  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]

  ##### Weighting object for diagnostics #####
  weight_object = weighting(outcome, treatment, trial, selection_covariates, data, selection_method, is_data_disjoint, trim_pop)

  participation_probs = weight_object$participation_probs
  weights = weight_object$weights
  g_index = gen_index(participation_probs$probs_population, participation_probs$probs_trial)

  ##### Generalize results #####
  ## First, estimate SATE
=======

  ##### just keep the data we need #####
  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]

  ##### trim population #####
  if(trim_pop == TRUE){
    data = trim_pop(trial, covariates, data)
  }

  ##### estimate SATE #####
>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
  SATE_model = lm(as.formula(paste(outcome,treatment,sep="~")), data = data)

  SATE = summary(SATE_model)$coefficients[treatment, "Estimate"]
  SATE_se = summary(SATE_model)$coefficients[treatment, "Std. Error"]

  SATE_CI_l = SATE - 1.96*SATE_se
  SATE_CI_u = SATE + 1.96*SATE_se

  SATE_results = c(SATE,SATE_se,SATE_CI_l,SATE_CI_u)

<<<<<<< HEAD
  ## Weighting results
  if(method == "weighting"){
    TATE_results = weight_object$TATE
=======
  ##### Weighting Methods #####
  if(method == "weighting"){
    TATE_results = weighting(outcome, treatment, trial, selection_covariates, data, method = selection_method, is_data_disjoint = is_data_disjoint)
>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
  }

  ## BART results
  if(method == "BART"){
    TATE_results = "NOT READY YET"
  }

  ## TMLE results
  if(method == "TMLE"){
<<<<<<< HEAD
    TATE_results = tmle(outcome, treatment, trial, selection_covariates, data)$TATE
  }

  ## put together results table
=======
    TATE_results = tmle(outcome, treatment, trial, selection_covariates, data)
  }

>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
  result.tab = rbind(SATE_results, TATE_results)
  colnames(result.tab) = c("Estimate","Std. Error","95% CI Lower","95% CI Upper")
  row.names(result.tab) = c("SATE","TATE")

<<<<<<< HEAD
  ##### Items to save to "generalize" object #####
  out = list(
    result.tab = result.tab,
    n_excluded = n_excluded,
    g_index = g_index,
    outcome = outcome,
    treatment = treatment,
    trial = trial,
    selection_covariates = selection_covariates,
    method = method,
    selection_method = selection_method,
    trim_pop = trim_pop
  )

  class(out) = "generalize"

  return(out)
  #invisible(out) --> returns output but doesn't print anything
=======
  return(result.tab)
>>>>>>> 9556e35d29cc067f2a5261ef6c3024119ac2d1d4
}

# S3 Methods (roxygenize voodoo)
# summary.generalize #first parameter needs to be called "object"
# print.generalize #first parameter needs to be called "x" (maybe use "invisible(x)")
# print.summary.generalize
