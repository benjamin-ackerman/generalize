# outcome: name of outcome variable
# treatment: name of treatment variable
# selection_formula: probability of trial participation formula
#

# method: c("weighting","bart","tmle")
# weighting_method: c("lr","rf","lasso","sl")

##### package dependencies: tmle


generalize <- function(outcome, treatment, selection_formula, data, method, weight_method = "lr",outcome_formula = NULL){

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
  data= data[rownames(na.omit(data[,all.vars(selection_formula)])),c(outcome, treatment, all.vars(selection_formula))]

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

    data = trim_pop(selection_formula, data)
    selection.model <- glm(selection_formula,family="quasibinomial",data=data)

    data_new0 <- data
    data_new1 <- data
    data_new0$a <- 0
    data_new1$a <- 1

    data$a1 <- predict(selection.model, newdata=data_new1, type="response")
    data$a0 <- predict(selection.model, newdata=data_new0, type="response")

    #arbitrarily assign A to those with missing A (those not in trial)
    tmle.model0 <- tmle::tmle(Y=data[,outcome],
                        A=ifelse(!is.na(data[,treatment]),data[,treatment],0),
                        W=data[,covariates],
                        Delta=data[,trial_membership],
                        g1W=0.5,
                        family='binomial',
                        gbound=c(0,1),
                        pDelta1=cbind(data$a0,data$a1)
    )
    OR = round(tmle.model0$estimate$OR$psi,2)
    CI = paste0("(",round(tmle.model0$estimate$OR$CI[1],2),
                "-",round(tmle.model0$estimate$OR$CI[2],2),")")

    results = list(TATE = OR, TATE_CI = CI)
  }

  return(results)
}
