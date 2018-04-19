#' BART for estimating TATE
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.
#' @param seed numeric. By default, the seed is set to 13783, otherwise can be specified (such as for simulation purposes).
#' @return \code{generalize_bart} returns a list of the TATE estimate, standard error, and 95\% CI bounds

generalize_bart <- function(outcome, treatment, trial, selection_covariates, data,is_data_disjoint = TRUE,seed){

  ##### set the seed #####
  if(missing(seed)){
    seed = 13783
  }

  set.seed(seed)

  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]

  # Prepare target population test set
  xtest = data[which(data[,trial] == 0),selection_covariates]

  ##### if the data are not disjoint:
  if(is_data_disjoint == FALSE){
    xtest = data[,selection_covariates]
  }

  xtest = rbind(xtest,xtest)
  xtest[,treatment] = rep(c(1,0),each = nrow(xtest)/2)

  # Training data based on trial
  xtrain = data[which(data[,trial] == 1 & !is.na(data[,outcome])),c(treatment,selection_covariates)]
  ytrain = data[which(data[,trial] == 1 & !is.na(data[,outcome])),outcome]

  bart.out <- BayesTree::bart(x.train = xtrain,
                              y.train = ytrain,
                              keeptrainfits = FALSE,
                              x.test = xtest,
                              verbose = FALSE)

  if(dim(table(data[,outcome])) == 2){
    bart.fits <- pnorm(bart.out$yhat.test[,1:(nrow(xtest)/2)]) -
      pnorm(bart.out$yhat.test[,(1+(nrow(xtest)/2)):nrow(xtest)])
  }
  else{
    bart.fits <- bart.out$yhat.test[,1:(nrow(xtest)/2)] -
      bart.out$yhat.test[,(1+(nrow(xtest)/2)):nrow(xtest)]
  }

  TATE = mean(apply(bart.fits,1,mean))
  TATE_se = sd(apply(bart.fits,1,mean))
  TATE_CI_l = TATE - 1.96*TATE_se
  TATE_CI_u = TATE + 1.96*TATE_se

  # #Function to get ORs from posterior draws after running BART
  # getORs = function(i){
  #
  #   pred_probs=pnorm(bart.out$yhat.test[i,])
  #   #pred_outcomes = sapply(pred_probs,function(x) rbinom(1,1,x))
  #
  #   placebo=pred_probs[which(xtest$treat == 0)]
  #   trt=pred_probs[which(!xtest$treat == 0)]
  #
  #   placebo_prob = mean(placebo,na.rm=TRUE)
  #   treatment_prob = mean(trt,na.rm=TRUE)
  #
  #   (treatment_prob/(1-treatment_prob))/(placebo_prob/(1-placebo_prob))
  #
  #   #(sum(trt)/(length(trt)-sum(trt)))/(sum(placebo)/(length(placebo)-sum(placebo)))
  # }

  # ORs = sapply(1:ndpost,getORs)

  TATE = list(estimate = TATE, se = TATE_se, CI_l = TATE_CI_l, CI_u = TATE_CI_u)

  out = list(TATE = TATE)

  return(out)
}
