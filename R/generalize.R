### MAKE BART, TMLE their own function --> makes things very modular, can add new methods
### separate out sample indicator and covariates? instead of formula?
### get rid of outcome model parameter
### instead of selection_formula, separate selection indicator and covariates

#' Estimate weights for generalizing ATE by predicting probability of trial participation
#'
#' @param outcome variable denoting outcome
#' @param treatment variable denoting binary treatment assignment (ok if only available in trial, not population)
#' @param selection_formula an object of class "formula." The formula specifying the model for trial participation.  Lefthand side should be a binary variable indicating trial membership, and righthand side should contain pre-treatment covariates measured in data set.
#' @param data a data frame containing the variables specified in the model
#' @param method choose method to generalize average treatment effect.  Default is "weighting" (weighting by the odds of participation probability).  Other methods supported are "bart" (Bayesian Additive Regression Trees - NOT READY YET) and "tmle" (Targeted Maximum Likelihood Estimation)
#' @param weight_method choose method to predict the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are random forests ("rf") and lasso ("lasso")
#' @param outcome_formula an object of class "formula." Can specify an optional outcome model to include pre-treatment covariates.
#' @return \code{generalize} returns an object of the class "generalize", containing the following: \code{TATE} (target population average treatment effect), \code{TATE_CI} (95% Confidence Interval for TATE).  If outcome is binary, reports TATE as risk difference as well as odds ratio, with accompanying CIs
#' @examples
#' generalize(outcome = "STUDYCOMPLETE", treatment = "treat", selection_formula = trial ~ age + sex + race, data = ctn_data, method = "weighting")
#' generalize(outcome = "STUDYCOMPLETE", treatment = "treat", selection_formula = trial ~ age + sex + race, data = ctn_data, method = "tmle")


generalize <- function(outcome, treatment, selection_formula, data,
                       method = "weighting", weight_method = "lr",outcome_formula = NULL){

  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  if(anyNA(match(outcome,names(data)))){
    stop("Outcome is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(treatment,names(data)))){
    stop("Treatment is not a variable in the data provided!",call. = FALSE)
  }

  if(class(selection_formula) != "formula" | (!is.null(outcome_formula) & class(outcome_formula) != "formula")){
    stop("Must enter a valid formula!",call. = FALSE)
  }


  if(anyNA(match(all.vars(selection_formula),names(data)))){
    missing_variables = all.vars(selection_formula)[is.na(match(all.vars(selection_formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

  if(!is.null(outcome_formula) & anyNA(match(all.vars(outcome_formula),names(data)))){
    missing_variables = all.vars(outcome_formula)[is.na(match(all.vars(outcome_formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

  trial_membership = all.vars(selection_formula)[1]
  covariates = all.vars(selection_formula)[-1]

  if(!length(unique(data[,trial_membership])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(!method %in% c("weighting","BART","TMLE")){
    stop("Invalid method!",call. = FALSE)
  }

  ##### just keep the data we need #####
  data = data[rownames(na.omit(data[,all.vars(selection_formula)])),c(outcome, treatment, all.vars(selection_formula))]

  ##### Weighting Methods #####
  if(method == "weighting"){
    data$weights = gen_weights(selection_formula, data = data, method = weight_method)$weights

    if(is.null(outcome_formula)){
      TATE_model = lm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights)
    }

    if(!is.null(outcome_formula)){
      TATE_model = lm(outcome_formula, data = data, weights = weights)
    }

    TATE = summary(TATE_model)$coefficients[treatment,"Estimate"]
    TATE_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]

    results = list(TATE = TATE, TATE_CI = c(TATE - 1.96*TATE_se, TATE + 1.96*TATE_se))

    if(dim(table(data[,outcome])) == 2){
      if(is.null(outcome_formula)){
        TATE_OR_model = glm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights,family='quasibinomial')
      }

      if(!is.null(outcome_formula)){
        TATE_OR_model = glm(outcome_formula, data = data, weights = weights,family='quasibinomial')
      }

      TATE_logOR = summary(TATE_model)$coefficients[treatment,"Estimate"]
      TATE_logOR_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]

      results = list(TATE = TATE, TATE_CI = c(TATE - 1.96*TATE_se, TATE + 1.96*TATE_se), TATE_OR = exp(TATE_logOR), TATE_OR_CI = c(exp(TATE_logOR - 1.96*TATE_logOR_se),exp(TATE_logOR + 1.96*TATE_logOR_se)))
    }
  }

  if(method == "BART"){
    results = "NOT READY YET"
  }

  if(method == "TMLE"){
    results = tmle(outcome, treatment, selection_formula, data)
  }

  return(results)
}
