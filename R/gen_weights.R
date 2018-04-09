## DONE: Parameter to specify what kind of weight you want (specifying to the "S = 1" or "S = 0 and 1" group) not disjoint: weight by inverse probability
## Separate the two tasks: 1 is create the weights, 1 is assessing similarities
## Don't print anything from this, have it be the behind the scenes
## Can use this to send to diagnostics, G-index, or generalize function
## Diagnostics: covariate balance, weighting --> include some sort of density plot

#' Estimate weights for generalizing ATE by predicting probability of trial participation
#'
#' @param formula an object of class "formula". The formula specifying the model for trial participation.  Lefthand side should be a binary variable indicating trial membership, and righthand side should contain pre-treatment covariates measured in data set.
#' @param data a data frame containing the variables specified in the model
#' @param method choose method to predict the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are random forests ("rf") and lasso ("lasso")
#' @return \code{gen_weights} returns an object of the class "gen_weights", containing the following: \code{participation_probs} (predicted probabilities of trial participation), \code{weights} (weights constructed from predicted probabilities - see description for more details), \code{gen_index} (generalizability index)
#' @examples
#' gen_weights(trial ~ age + sex + race, data = ctn_data)
#' gen_weights(trial ~ age + sex + race, data = ctn_data, method = 'rf')

gen_weights = function(formula, data, method = "lr", is.data.disjoint = TRUE){

  ### Get variable names from the formula ###
  trial_membership = all.vars(formula)[1]
  covariates = all.vars(formula)[-1]

  ### Checks ###

  # if (missing(data)) {
  #   data <- environment(formula)
  #   #stop("Data must be specified.", call. = FALSE)
  # }

  if(class(formula) != "formula"){
    stop("Must enter a valid formula!",call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if(anyNA(match(all.vars(formula),names(data)))){
    missing_variables = all.vars(formula)[is.na(match(all.vars(formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," not in the data provided"))
  }

  if(!length(unique(data[,trial_membership])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }


  ### Clean up data from missing values ###
  data = data[rownames(na.omit(data[,all.vars(formula)])),]

  ### Logistic Regression

  if(method == "lr"){
    ps = predict(glm(formula, data = data, family='quasibinomial'),type = 'response')
  }

  if(method == "rf"){
    formula = as.formula(paste0("as.factor(",trial_membership,") ~ ",as.character(formula)[3]))
    ps = predict(randomForest::randomForest(formula, data=data, na.action=na.omit, sampsize = 454, ntree=1500),type = 'prob')[,2]
  }

  if(method == "lasso"){
    test.x = model.matrix(~ -1 + ., data=data[,covariates])
    test.y = data[,trial_membership]
    ps = as.numeric(predict(glmnet::cv.glmnet(
      x=test.x,
      y=test.y,
      family="binomial"
    ),newx=test.x,s="lambda.1se",type="response"))
  }

  if(is.data.disjoint == TRUE){
    weights = ifelse(data[,trial_membership]==0,0,ps/(1-ps))
  }

  if(is.data.disjoint == FALSE){
    weights = ifelse(data[,trial_membership]==0,0,1/ps)
  }

  participation_probs = list(probs_population = ps[which(data[,trial_membership]==0)],
                     probs_trial = ps[which(data[,trial_membership]==0)])

  out = list(
    participation_probs = participation_probs,
    weights = weights
  )

  return(out)
}


