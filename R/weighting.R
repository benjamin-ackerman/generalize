#' Weighting method
#'
#' @param outcome variable denoting outcome
#' @param treatment variable denoting binary treatment assignment (ok if only available in trial, not population)
#' @param selection_formula an object of class "formula." The formula specifying the model for trial participation.  Lefthand side should be a binary variable indicating trial membership, and righthand side should contain pre-treatment covariates measured in data set.
#' @param data a data frame containing the variables specified in the model
#' @param selection_method choose method to predict the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are random forests ("rf") and lasso ("lasso")
#' @param outcome_formula an object of class "formula." Can specify an optional outcome model to include pre-treatment covariates.
#' @return \code{generalize} returns an object of the class "generalize", containing the following: \code{TATE} (target population average treatment effect), \code{TATE_CI} (95% Confidence Interval for TATE).  If outcome is binary, reports TATE as risk difference as well as odds ratio, with accompanying CIs
#' @examples
#' generalize(outcome = "STUDYCOMPLETE", treatment = "treat", selection_formula = trial ~ age + sex + race, data = ctn_data, method = "weighting")
#' generalize(outcome = "STUDYCOMPLETE", treatment = "treat", selection_formula = trial ~ age + sex + race, data = ctn_data, method = "tmle")

weighting <- function(outcome, # variable name
                       treatment, # variable name: must be binary indicator of treatment, values can be missing in population data
                       trial, # variable name: must be binary indicator
                       selection_covariates, # a vector of covariate names in data set that predict trial participation
                       data, # data frame containing data
                       method = "lr", # weighting method, can either be "lr" (logistic regression), "rf" (random forests), "lasso" (lasso)
                       is.data.disjoint = TRUE # logical to determine how to calculate weights
                       ){

  ##### just keep the data we need #####
  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]

  data$weights = gen_weights(trial, selection_covariates, data, method = selection_method, is.data.disjoint)$weights

  TATE_model = lm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights)

  TATE = summary(TATE_model)$coefficients[treatment,"Estimate"]
  TATE_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]

  TATE_CI_l = TATE - 1.96*TATE_se
  TATE_CI_u = TATE + 1.96*TATE_se

  out = c(TATE, TATE_se, TATE_CI_l, TATE_CI_u)

  ##### DEALING WITH BINARY OUTCOMES #####
  # if(dim(table(data[,outcome])) == 2){
  #   if(is.null(outcome_formula)){
  #     TATE_OR_model = glm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights,family='quasibinomial')
  #   }
  #
  #   if(!is.null(outcome_formula)){
  #     TATE_OR_model = glm(outcome_formula, data = data, weights = weights,family='quasibinomial')
  #   }
  #
  #   TATE_logOR = summary(TATE_model)$coefficients[treatment,"Estimate"]
  #   TATE_logOR_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]
  #
  #   results = list(TATE = TATE,
  #                  TATE_CI = c(TATE - 1.96*TATE_se, TATE + 1.96*TATE_se),
  #                  TATE_OR = exp(TATE_logOR),
  #                  TATE_OR_CI = c(exp(TATE_logOR - 1.96*TATE_logOR_se),exp(TATE_logOR + 1.96*TATE_logOR_se)))
  # }

  return(out)
}
