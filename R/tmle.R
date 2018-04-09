tmle <- function(outcome,treatment, formula, data){
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  if(anyNA(match(outcome,names(data)))){
    stop("Outcome is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(treatment,names(data)))){
    stop("Treatment is not a variable in the data provided!",call. = FALSE)
  }

  if(class(formula) != "formula"){
    stop("Must enter a valid formula!",call. = FALSE)
  }

  if(anyNA(match(all.vars(selection_formula),names(data)))){
    missing_variables = all.vars(selection_formula)[is.na(match(all.vars(selection_formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

  trial_membership = all.vars(selection_formula)[1]
  covariates = all.vars(selection_formula)[-1]

  data = trim_pop(formula, data)$trimmed_data
  selection.model <- glm(formula,family="quasibinomial",data=data)

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
