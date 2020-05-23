#' Estimate weights for generalizing ATE by predicting probability of trial participation
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf"), Lasso ("lasso"), GBM ("gbm") and SuperLearner ("super").
#' @param sl_library vector of SuperLearner library methods. If `selection_method` = 'super', specify names of methods to include in library. Default is NULL.
#' @param survey_weights variable name of population data's complex survey weights. Default is FALSE: if FALSE, then population data do not come a complex survey and weights do not need to be incorporated in estimation. NOTE: SURVEY WEIGHTS ONLY SUPPORTED FOR USE WHEN USING WEIGHTING METHODS FOR THE TIME BEING. Survey weights will not be incorporated in BART or TMLE.
#' @param trim_weights logical. If TRUE, then trim the weights to the value specified in `trim_pctile`. Default is FALSE.
#' @param trim_pctile numeric. If `trim_weights` is TRUE, then specify what percentile weights should be trimmed to. Default is 0.97.
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights - see details for more information.
#' @param seed numeric. By default, the seed is set to 13783, otherwise can be specified (such as for simulation purposes).

#' @export
weighting = function(outcome, treatment, trial, selection_covariates, data,
                     selection_method = "lr", sl_library = NULL, survey_weights = FALSE, trim_weights=FALSE, trim_pctile = .97, is_data_disjoint = TRUE,seed){
  . = n_persons  = TATE_se = NULL
  rm(list = c( "n_persons", "TATE_se", "."))

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

  if(!selection_method %in% c("lr","rf","lasso","gbm","super")){
    stop("Invalid method!",call. = FALSE)
  }

  ### Clean up data from missing values ###

  if(survey_weights == FALSE){
    data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, selection_covariates)]
    data$s_weights = 1
  } else{
    data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),c(outcome, treatment, trial, survey_weights, selection_covariates)]
    data$s_weights = ifelse(data[,trial] == 1, 1, data[,survey_weights])

    ### NORMALIZE THE SURVEY WEIGHTS TO SUM TO THE NUMBER OF SURVEY PARTICIPANTS:
    normalize_factor = mean(data$s_weights[which(data[,trial] == 0)], na.rm = TRUE)
    data$s_weights[which(data[,trial] == 0)] = data$s_weights[which(data[,trial] == 0)]/normalize_factor

    if(selection_method == "rf"){
      data = data %>%
        dplyr::filter(get(trial) == 0) %>%
        dplyr::mutate(n_persons = ceiling(get(survey_weights))) %>%
        tidyr::uncount(n_persons) %>%
        dplyr::bind_rows(data %>% filter(get(trial) == 1))
    }
  }

  # Change selection method to fit WeightIt

  if(selection_method %in% c("lr","gbm","super")){
    if(selection_method == "lr"){selection_method = "ps"}

    formula = as.formula(paste(trial, paste(selection_covariates,collapse="+"),sep="~"))
    ps = WeightIt::weightit(formula, data = data, method = selection_method,
                  estimand="ATT",focal="0",s.weights = data$s_weights,
                  stop.method="ks.mean",#gbm parameters
                  SL.library = sl_library)$ps #superlearner parameters
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
      weights = data[,"s_weights"],
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
  if(trim_weights){
    cutoff = as.numeric(quantile(data$weights[which(data[,trial] == 1)], trim_pctile))
    data$weights[which(data[,trial] == 1)] = ifelse(data$weights[which(data[,trial] == 1)] > cutoff, cutoff, data$weights[which(data[,trial] == 1)])
  }

  participation_probs = list(population = ps[which(data[,trial]==0)],
                             trial = ps[which(data[,trial]==1)])

  if(is.null(outcome) & is.null(treatment)){TATE = NULL}
  else{
    ##### ESTIMATE POPULATION AVERAGE TREATMENT EFFECT #####
    formula = as.formula(paste(outcome,treatment,sep="~"))

    ATE_design = survey::svydesign(id = ~1, data = data[which(data[,trial] == 1),], weights = data$weights[which(data[,trial] == 1)])

    if(length(table(data[,outcome]))!=2){
      model = survey::svyglm(formula, design = ATE_design, family='gaussian')
      TATE = summary(model)$coefficients[treatment,"Estimate"]
      TATE_se = summary(model)$coefficients[treatment,"Std. Error"]
      TATE_CI_l = as.numeric(confint(model)[treatment,])[1]
      TATE_CI_u =as.numeric(confint(model)[treatment,])[2]
    } else{
      model = survey::svyglm(formula, design = ATE_design, family='quasibinomial')
      TATE = exp(summary(model)$coefficients[treatment,"Estimate"])
      TATE_se = NA
      TATE_CI_l = as.numeric(exp(confint(model)[treatment,]))[1]
      TATE_CI_u = as.numeric(exp(confint(model)[treatment,]))[2]
    }

    TATE = list(estimate = TATE, se = TATE_se, CI_l = TATE_CI_l, CI_u = TATE_CI_u)
  }

  ##### Items to return out #####
  out = list(participation_probs = participation_probs,
             weights = data$weights,
             TATE = TATE)

  return(out)
}
