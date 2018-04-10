## DONE: Parameter to specify what kind of weight you want (specifying to the "S = 1" or "S = 0 and 1" group) not disjoint: weight by inverse probability
## Separate the two tasks: 1 is create the weights, 1 is assessing similarities
## Don't print anything from this, have it be the behind the scenes
## Can use this to send to diagnostics, G-index, or generalize function
## Diagnostics: covariate balance, weighting --> include some sort of density plot

#' Estimate weights for generalizing ATE by predicting probability of trial participation
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param trim_pop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of trial covariates.
#' @return \code{gen_weights} returns an object of the class "gen_weights", containing the following: \code{participation_probs} (predicted probabilities of trial participation), \code{weights} (weights constructed from predicted probabilities - see description for more details), \code{gen_index} (generalizability index)
#' @examples
#' gen_weights(trial ~ age + sex + race, data = ctn_data)
#' gen_weights(trial ~ age + sex + race, data = ctn_data, method = 'rf')

gen_weights = function(trial, selection_covariates, data, selection_method = "lr",
                       is_data_disjoint = TRUE, trim_pop = TRUE){

  ### Make input method lower case ###
  selection_method = tolower(selection_method)

  ### Checks ###
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if(anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(unique(data[,trial])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid method!",call. = FALSE)
  }

  ### Clean up data from missing values ###
  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),]

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

  ### Generate Weights ###
  if(is_data_disjoint == TRUE){
    weights = ifelse(data[,trial]==0,0,ps/(1-ps))
  }

  if(is_data_disjoint == FALSE){
    weights = ifelse(data[,trial]==0,0,1/ps)
  }

  participation_probs = list(probs_population = ps[which(data[,trial]==0)],
                             probs_trial = ps[which(data[,trial]==0)])

  out = list(method = selection_method,
             participation_probs = participation_probs,
             weights = weights)

  return(out)
}


