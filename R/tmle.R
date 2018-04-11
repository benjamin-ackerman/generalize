#' TMLE
#'
#' @param outcome variable denoting outcome
#' @param treatment variable denoting binary treatment assignment (ok if only available in trial, not population)
#' @param data a data frame containing the variables specified in the model
#' @param method choose method to generalize average treatment effect.  Default is "weighting" (weighting by the odds of participation probability).  Other methods supported are "bart" (Bayesian Additive Regression Trees - NOT READY YET) and "tmle" (Targeted Maximum Likelihood Estimation)
#' @param weight_method choose method to predict the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are random forests ("rf") and lasso ("lasso")
#' @param outcome_formula an object of class "formula." Can specify an optional outcome model to include pre-treatment covariates.
#' @return \code{generalize} returns an object of the class "generalize", containing the following: \code{TATE} (target population average treatment effect), \code{TATE_CI} (95% Confidence Interval for TATE).  If outcome is binary, reports TATE as risk difference as well as odds ratio, with accompanying CIs
#' @examples

tmle <- function(outcome, # variable name
                 treatment, # variable name: must be binary indicator of treatment, values can be missing in population data
                 trial, # variable name: must be binary indicator
                 selection_covariates, # a vector of covariate names in data set that predict trial participation
                 data # data frame containing data
                 ){

  data = trim_pop(trial = trial, covariates = selection_covariates, data = data)$trimmed_data

  formula = as.formula(paste(trial, paste(selection_covariates,collapse="+"),sep="~"))

  selection.model <- glm(formula,family="quasibinomial",data=data)

  data_new0 <- data
  data_new1 <- data
  data_new0$a <- 0
  data_new1$a <- 1

  data$a1 <- predict(selection.model, newdata=data_new1, type="response")
  data$a0 <- predict(selection.model, newdata=data_new0, type="response")

  tmle.family = "gaussian"

  ### Trying to account for binary outcomes
  # if(dim(table(data[,outcome])) == 2){
  #   tmle.family = "binomial"
  # }

  tmle.model <- tmle::tmle(Y=data[,outcome],
                              A=ifelse(!is.na(data[,treatment]),data[,treatment],0),
                              W=data[,covariates],
                              Delta=data[,trial_membership],
                              g1W=0.5,
                              family= tmle.family,
                              gbound=c(0,1),
                              pDelta1=cbind(data$a0,data$a1)
                             )

  TATE = tmle.model$estimates$ATE$psi
  TATE_se = sqrt(tmle.model$estimates$ATE$var.psi)
  TATE_CI_l = TATE - 1.96*TATE_se
  TATE_CI_u = TATE + 1.96*TATE_se

  # TATE_OR = round(tmle.model$estimate$OR$psi,2)
  # TATE_se

  # CI = paste0("(",round(tmle.model$estimate$OR$CI[1],2),
  #             "-",round(tmle.model$estimate$OR$CI[2],2),")")

  out = list(TATE = c(TATE, TATE_se, TATE_CI_l, TATE_CI_u))

  return(out)
}
