#' Estimate weights for generalizing ATE by predicting probability of trial participation
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param seed numeric. By default, the seed is set to 13783, otherwise can be specified (such as for simulation purposes).

weighting = function(outcome, treatment, trial, selection_covariates, data,
                     selection_method = "lr", is_data_disjoint = TRUE,seed){

  ##### set the seed #####
  if(missing(seed)){
    seed = 13783
  }

  set.seed(seed)

  ### Make input method lower case ###
  selection_method = tolower(selection_method)

  ### Checks ###
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if(anyNA(match(outcome,names(data)))){
    stop("Outcome is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(treatment,names(data)))){
    stop("Treatment is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(na.omit(unique(data[,trial]))) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid method!",call. = FALSE)
  }

  ### Clean up data from missing values ###
  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]

  ### Generate Participation Probabilities ###
  # Logistic Regression
  if(selection_method == "lr"){
    formula = as.formula(paste(trial, paste(selection_covariates,collapse="+"),sep="~"))
    ps = predict(glm(formula, data = data, family='quasibinomial'),type = 'response')
  }

  # Random Forests
  if(selection_method == "rf"){
    formula = as.formula(paste( paste("as.factor(",trial,")"), paste(selection_covariates,collapse="+"),sep="~"))
    ps = predict(randomForest::randomForest(formula, data=data, na.action=na.omit, sampsize = 454, ntree=1500),type = 'prob')[,2]
  }

  # Lasso
  if(selection_method == "lasso"){
    test.x = model.matrix(~ -1 + ., data=data[,selection_covariates])
    test.y = data[,trial]
    ps = as.numeric(predict(glmnet::cv.glmnet(
      x=test.x,
      y=test.y,
      family="binomial"
    ),newx=test.x,s="lambda.1se",type="response"))
  }

  ### Set any participation probabilities of 0 in the trial to the minimum non-zero value ###
  if(any(ps[which(data[,trial]==1)] == 0)){
    ps[which(data[,trial] == 1 & ps == 0)] = min(ps[which(data[,trial] == 1 & ps != 0)], na.rm=TRUE)
  }

  ### Generate Weights ###
  if(is_data_disjoint == TRUE){
    data$weights = ifelse(data[,trial]==0,0,(1-ps)/ps)
  }

  if(is_data_disjoint == FALSE){
    data$weights = ifelse(data[,trial]==0,0,1/ps)
  }

  # Trim any of the weights if necessary
  data$weights[which(data$weights == 0 & data[,trial] == 1)] = quantile(data$weights[which(data[,trial]==1)],0.01,na.rm=TRUE)

  participation_probs = list(population = ps[which(data[,trial]==0)],
                             trial = ps[which(data[,trial]==1)])

  if(is.null(outcome) & is.null(treatment)){TATE = NULL}
  else{

  ##### ESTIMATE POPULATION AVERAGE TREATMENT EFFECT #####
  TATE_model = lm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights)

  TATE = summary(TATE_model)$coefficients[treatment,"Estimate"]
  TATE_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]

  TATE_CI_l = TATE - 1.96*TATE_se
  TATE_CI_u = TATE + 1.96*TATE_se

  TATE = list(estimate = TATE, se = TATE_se, CI_l = TATE_CI_l, CI_u = TATE_CI_u)
  }

  ##### Items to return out #####
  out = list(participation_probs = participation_probs,
             weights = data$weights,
             TATE = TATE)

  return(out)
}
