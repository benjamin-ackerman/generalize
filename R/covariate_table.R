#' Create Covariate Balance Table
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data

covariate_table <- function(trial, selection_covariates, data, weighted_table = FALSE, selection_method = "lr", is_data_disjoint = TRUE){

  if(weighted_table == FALSE){
    data = data %>%
      tidyr::drop_na(selection_covariates) %>%
      as.data.frame()

    dmy <- caret::dummyVars(" ~ .", data = data[,selection_covariates], sep="_")
    expanded.data = data.frame(trial = data[, trial], predict(dmy, newdata = data[,selection_covariates]))
    names(expanded.data)=stringr::str_replace_all(names(expanded.data),"\\.","-")

    means.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1], mean) %>%
      t() %>% as.data.frame()

    names(means.tab)[which(means.tab["trial",]=="1")] = "trial"
    names(means.tab)[which(means.tab["trial",]=="0")] = "population"

    means.tab = means.tab[!rownames(means.tab) %in% c("trial","X.Intercept."),]

    n_trial = as.numeric(table(expanded.data$trial)["1"])
    n_pop = as.numeric(table(expanded.data$trial)["0"])

    sd.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_all(var) %>% t() %>% as.data.frame()

    sd.tab$pooled_sd = sqrt(((n_trial - 1)*sd.tab$V1 + (n_pop - 1)*sd.tab$V2)/(n_trial + n_pop - 2))

    names(sd.tab)[which(sd.tab["trial",]=="1")] = "trial_var"
    names(sd.tab)[which(sd.tab["trial",]=="0")] = "population_var"
    names(sd.tab)[3] = "pooled_sd"

    sd.tab = sd.tab[!rownames(sd.tab) %in% c("trial", "X.Intercept."),]

    names(sd.tab) = c("trial_var", "population_var", "pooled_sd")
  }

  if(weighted_table == TRUE){
    data = data %>%
      tidyr::drop_na(selection_covariates) %>%
      as.data.frame()

    data$weights = weighting(outcome = NULL, treatment = NULL, trial = trial,
                             selection_covariates = selection_covariates, data = data,
                             selection_method = selection_method, is_data_disjoint = is_data_disjoint)$weights
    data$weights = ifelse(data[,trial] == 0, 1, data$weights)

    dmy <- caret::dummyVars(" ~ .", data = data[,selection_covariates], sep="_")
    expanded.data = data.frame(trial = data[, c(trial)], predict(dmy, newdata = data[,selection_covariates]), weights = data$weights)
    names(expanded.data)=stringr::str_replace_all(names(expanded.data),"\\.","-")

    means.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1], funs(weighted.mean(., weights))) %>%
      dplyr::select(-`weights`) %>%
      t() %>%
      as.data.frame()

    names(means.tab)[which(means.tab["trial",]=="1")] = "trial"
    names(means.tab)[which(means.tab["trial",]=="0")] = "population"

    means.tab = means.tab[!rownames(means.tab) %in% c("trial", "X.Intercept."),]

    n_trial = as.numeric(table(expanded.data$trial)["1"])
    n_pop = as.numeric(table(expanded.data$trial)["0"])

    sd.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[which(!names(expanded.data)%in% c("trial","X.Intercept."))],
                          funs(sum(weights * (. - weighted.mean(.,weights))^2)/sum(weights))) %>%
      dplyr::select(-`weights`) %>%
      t() %>%
      as.data.frame()

    sd.tab$pooled_sd = sqrt(((n_trial - 1)*sd.tab$V1 + (n_pop - 1)*sd.tab$V2)/(n_trial + n_pop - 2))

    names(sd.tab)[which(sd.tab["trial",]=="1")] = "trial_var"
    names(sd.tab)[which(sd.tab["trial",]=="0")] = "population_var"
    names(sd.tab)[3] = "pooled_sd"

    sd.tab = sd.tab[!rownames(sd.tab) %in% c("trial", "X.Intercept."),]
  }

  covariate_table = means.tab %>%
    dplyr::bind_cols(sd.tab) %>%
    dplyr::mutate(ASMD = round(abs((trial - population)/pooled_sd),3)) %>%
    dplyr::select(trial, population, ASMD)

  if(weighted_table == FALSE){
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial","X.Intercept."))
  }

  if(weighted_table == TRUE){
    names(covariate_table)[which(names(covariate_table) == "trial")] = "trial (weighted)"
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial","X.Intercept.","weights"))
  }

  return(covariate_table)
}
