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

    expanded.data = data.frame(trial = data[, trial],
                               model.matrix(~-1 + ., data = data[, selection_covariates]))

    means.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1], mean) %>%
      t() %>% as.data.frame()

    means.tab = means.tab[-1,]

    names(means.tab) = c("trial", "population")
    n_trial = as.numeric(table(expanded.data[,"trial"]))[2]
    n_pop = as.numeric(table(expanded.data[,"trial"]))[1]
    sd.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_all(var) %>% t() %>% as.data.frame() %>% .[-1,] %>%
      mutate(pooled_sd = sqrt(((n_trial - 1) * V1 + (n_pop - 1) * V2)/(n_trial + n_pop - 2)))
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

    expanded.data = data.frame(trial = data[,trial], model.matrix(~ -1 + ., data = data[,c(selection_covariates,"weights")]))

    means.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1], funs(weighted.mean(., weights))) %>%
      dplyr::select(-`weights`) %>%
      t() %>%
      as.data.frame() %>% .[-1,]

    names(means.tab) = c("trial","population")

    n_trial = as.numeric(table(expanded.data$trial))[2]
    n_pop = as.numeric(table(expanded.data$trial))[1]

    sd.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1],
                          funs(sum(weights * (. - weighted.mean(.,weights))^2)/sum(weights))) %>%
      dplyr::select(-`weights`) %>%
      t() %>%
      as.data.frame() %>% .[-1,] %>%
      mutate(pooled_sd = sqrt(((n_trial - 1)*V1 + (n_pop - 1)*V2)/(n_trial + n_pop - 2)))

    names(sd.tab) = c("trial_var","population_var","pooled_sd")
  }

  covariate_table = means.tab %>%
    dplyr::bind_cols(sd.tab) %>%
    dplyr::mutate(ASMD = round(abs((trial - population)/pooled_sd),3)) %>%
    dplyr::select(trial, population, ASMD)

  if(weighted_table == FALSE){
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial"))
  }

  if(weighted_table == TRUE){
    names(covariate_table)[1] = "trial (weighted)"
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial","weights"))
  }

  return(covariate_table)
}
