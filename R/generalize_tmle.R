#' TMLE for estimating TATE
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.
#' @param seed numeric. By default, the seed is set to 13783, otherwise can be specified (such as for simulation purposes).
#' @return \code{generalize_tmle} returns a list of the TATE estimate, standard error, and 95\% CI bounds

#' @export
generalize_tmle <- function(outcome, treatment, trial, selection_covariates, data,is_data_disjoint = TRUE,seed){

  ##### set the seed #####
  if(missing(seed)){
    seed = 13783
  }
  set.seed(seed)

  data = trim_pop(trial = trial, selection_covariates = selection_covariates, data = data)$trimmed_data

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
                           W=data[,selection_covariates],
                           Z=NULL,
                           Delta=data[,trial],
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

  TATE = list(estimate = TATE, se = TATE_se, CI_l = TATE_CI_l, CI_u = TATE_CI_u)

  out = list(TATE = TATE)

  return(out)
}
