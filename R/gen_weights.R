gen_weights = function(formula, data, method = "lr"){

  ### Get variable names from the formula ###
  trial_membership = all.vars(formula)[1]
  covariates = all.vars(formula)[-1]

  ### Checks ###

  # if (missing(data)) {
  #   data <- environment(formula)
  #   #stop("Data must be specified.", call. = FALSE)
  # }

  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if(anyNA(match(all.vars(formula),names(data)))){
    missing_variables = all.vars(formula)[is.na(match(all.vars(formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

  if(!length(unique(data[,trial_membership])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }


  ### Clean up data from missing values ###
  data = na.omit(data[,all.vars(formula)])

  ### Logistic Regression

  if(method == "lr"){
    ps = predict(glm(formula, data = data, family='quasibinomial'),type = 'response')
  }

  if(method == "rf"){
    formula = as.formula(paste0("as.factor(",trial_membership,") ~ ",as.character(formula)[3]))
    ps = predict(randomForest(formula, data=data, na.action=na.omit, sampsize = 454, ntree=1500),type = 'prob')[,2]
  }

  if(method == "lasso"){
    test.x = model.matrix(~ -1 + ., data=data %>% select(covariates))
    test.y = data[,trial_membership]
    ps = as.numeric(predict(cv.glmnet(
      x=test.x,
      y=test.y,
      family="binomial"
    ),newx=test.x,s="lambda.1se",type="response"))
  }

  weights = ifelse(data[,trial_membership]==0,0,ps/(1-ps))
  participation_probs = list(probs_population = ps[which(data[,trial_membership]==0)],
                     probs_trial = ps[which(data[,trial_membership]==0)])

  return(participation_probs)
}



gen_weights(trial ~ AGE + GENDER + RACE + HISPANIC, data = dat)

