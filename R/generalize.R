# Y, S, covariates, method

# Y: vector of outcomes
# S: vector of sample membership, binary
# X: matrix/data frame of pre-treatment covariates
# A: vector oftreatment assignment (can be missing in population)

# method: c("weighting","bart","tmle")
# weighting_method: c("lr","rf","lasso","sl")



generalize <- function(data, outcome, treatment, selection.formula, method = c("weighting","BART","TMLE"), weight_method = "lr",outcome.formula = NULL){

  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }


  if(anyNA(match(outcome,names(data)))){
    stop("Outcome is not a variable in the data provided!",call. = FALSE)
  }

  if(anyNA(match(treatment,names(data)))){
    stop("Treatment is not a variable in the data provided!",call. = FALSE)
  }

  if(class(selection.formula) != "formula" | (!is.null(outcome.formula) & class(outcome.formula) != "formula")){
    stop("Must enter a valid formula!",call. = FALSE)
  }


  if(anyNA(match(all.vars(selection.formula),names(data)))){
    missing_variables = all.vars(selection.formula)[is.na(match(all.vars(selection.formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

  if(!is.null(outcome.formula) & anyNA(match(all.vars(outcome.formula),names(data)))){
    missing_variables = all.vars(outcome.formula)[is.na(match(all.vars(outcome.formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

  trial_membership = all.vars(selection.formula)[1]
  covariates = all.vars(selection.formula)[-1]

  if(!length(unique(data[,trial_membership])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  ##### just keep the data we need #####
  data= data[rownames(na.omit(data[,all.vars(selection.formula)])),c(outcome, treatment, all.vars(selection.formula))]


  ##### generalize!!!

  ##### Weighting Methods #####
  if(method == "weighting"){
    data$weights = gen_weights(selection.formula, data = data, method = weight_method)$weights

    if(is.null(outcome.formula)){
      TATE_model = lm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights)
    }

    if(!is.null(outcome.formula)){
      TATE_model = lm(outcome.formula, data = data, weights = weights)
    }

    TATE = summary(TATE_model)$coefficients[treatment,"Estimate"]
    TATE_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]

    results = list(TATE = TATE, TATE_CI = c(TATE - 1.96*TATE_se, TATE + 1.96*TATE_se))

    if(dim(table(data[,outcome])) == 2){
      if(is.null(outcome.formula)){
        TATE_OR_model = glm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights,family='quasibinomial')
      }

      if(!is.null(outcome.formula)){
        TATE_OR_model = glm(outcome.formula, data = data, weights = weights,family='quasibinomial')
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
    results = "NOT READY YET"
  }

  return(results)
}


#########################################################################################################################


results = function(outcome){
  ### Unweighted Model (SATE) ###

  ### Weighted (GLM) Model (TATE) ###
  dat$weight = ifelse(dat$trial == 0, 0, dat$ps_lr/(1-dat$ps_lr))

  ## Outcome model ##
  weighted_glm_model = glm(get(outcome) ~ treat, family='quasibinomial',data = dat, weights = weight)

  OR = round(exp(summary(weighted_glm_model)$coefficients[-1,1]),2)
  CI = paste0("(",round(exp(confint.default(weighted_glm_model))[-1,1],2),                           "-",round(exp(confint.default(weighted_glm_model))[-1,2],2),")")
  weighted_glm = c(OR,CI)

  ### Weighted (RF) Model (TATE) ###
  dat$weight = ifelse(dat$trial == 0, 0, dat$ps_rf/(1-dat$ps_rf))

  ## Outcome model ##
  weighted_rf_model = glm(get(outcome) ~ treat, family='quasibinomial',data = dat, weights = weight)

  OR = round(exp(summary(weighted_rf_model)$coefficients[-1,1]),2)
  CI = paste0("(",round(exp(confint.default(weighted_rf_model))[-1,1],2),                           "-",round(exp(confint.default(weighted_rf_model))[-1,2],2),")")
  weighted_rf = c(OR,CI)

  ### Weighted (Lasso) Model (TATE) ###
  dat$weight = ifelse(dat$trial == 0, 0, dat$ps_lasso/(1-dat$ps_lasso))

  ## Outcome model ##
  weighted_lasso_model = glm(get(outcome) ~ treat, family='quasibinomial',data = dat, weights = weight)

  OR = round(exp(summary(weighted_lasso_model)$coefficients[-1,1]),2)
  CI = paste0("(",round(exp(confint.default(weighted_lasso_model))[-1,1],2),                           "-",round(exp(confint.default(weighted_lasso_model))[-1,2],2),")")
  weighted_lasso = c(OR,CI)

  ### Weighted (SuperLearner) Model (TATE) ###
  dat$weight = ifelse(dat$trial == 0, 0, dat$ps_sl/(1-dat$ps_sl))

  ## Outcome model ##
  weighted_sl_model = glm(get(outcome) ~ treat, family='quasibinomial',data = dat, weights = weight)

  OR = round(exp(summary(weighted_sl_model)$coefficients[-1,1]),2)
  CI = paste0("(",round(exp(confint.default(weighted_sl_model))[-1,1],2),                           "-",round(exp(confint.default(weighted_sl_model))[-1,2],2),")")
  weighted_sl = c(OR,CI)

  ### BART (TATE) ###
  bart_results = list.files()[str_detect(list.files(),paste0("BART_ORs",outcome))]

  bart_ORs = readRDS(bart_results)

  OR = round(mean(bart_ORs),2)
  CI = paste0("(",round(quantile(bart_ORs,prob=.025),2),"-",round(quantile(bart_ORs,prob=.975),2),")")
  bart = c(OR,CI)

  ### TMLE (TATE) ###
  tmle_dat = trim_data(dat)
  glm.model <- glm(trial ~ AGE + GENDER + RACE + HISPANIC + MARSTAT + EDUC + EMPLOY + METH7,family="quasibinomial",data=tmle_dat)

  data_new0 <- tmle_dat
  data_new1 <- tmle_dat
  data_new0$a <- 0
  data_new1$a <- 1

  tmle_dat$a1 <- predict(glm.model, newdata=data_new1, type="response")
  tmle_dat$a0 <- predict(glm.model, newdata=data_new0, type="response")

  #arbitrarily assign A to those with missing A (those not in trial)
  tmle.model0 <- tmle(Y=tmle_dat[,outcome],
                      A=ifelse(!is.na(tmle_dat$treat),tmle_dat$treat,0),
                      W=tmle_dat[,covariates],
                      Delta=tmle_dat$trial,
                      g1W=0.5,
                      family='binomial',
                      gbound=c(0,1),
                      pDelta1=cbind(tmle_dat$a0,tmle_dat$a1)
  )
  OR = round(tmle.model0$estimate$OR$psi,2)
  CI = paste0("(",round(tmle.model0$estimate$OR$CI[1],2),
              "-",round(tmle.model0$estimate$OR$CI[2],2),")")
  tmle = c(OR,CI)

  #### Return the outcomes ####
  result_tab = rbind(unweighted,weighted_glm,weighted_rf,weighted_lasso,weighted_sl,bart,tmle)
  colnames(result_tab) = c("OR","CI")
  rownames(result_tab) = NULL
  result_tab = as.data.frame(result_tab,stringsAsFactors=FALSE)
  result_tab$Method = c("Unweighted","Weighted (GLM)","Weighted (RF)","Weighted (Lasso)","Weighted (SL)","BART","TMLE")
  return(result_tab)
}
