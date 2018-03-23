# Y, S, covariates, method

# Y: vector of outcomes
# S: vector of sample membership, binary
# X: matrix/data frame of pre-treatment covariates
# A: vector oftreatment assignment (can be missing in population)

# method: c("weighting","bart","tmle")
# weighting_method: c("lr","rf","lasso","sl")



generalize <- function(Y, A, S, X, method){
  #data = data.frame(Y,A,S,X,stringsAsFactors = FALSE)

  if(tolower(method) == "weighting" & tolower(weighting_method) %in% c(NULL,"lr")){


  }


}


### read in: formula, data, method of weighting,
trial ~ x1 + x2 + x3 + ... + x4


formula = as.formula(study ~ x1 + x3 + x2*x4)
all.vars(formula)




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
