#' Create Covariate Balance Table
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data

covariate_table <- function(trial, selection_covariates, data){

  data = data[rownames(na.omit(data[,c(trial,selection_covariates)])),]

  expanded.data = data.frame(trial = data[,trial], model.matrix(~ 1 + ., data = data[,selection_covariates]))

  means.tab = expanded.data %>%
    dplyr::group_by(trial == 0) %>%
    dplyr::summarise_all(mean) %>%
    select(-trial,-`trial == 0`,-`X.Intercept.`) %>%
    t() %>%
    as.data.frame()

  names(means.tab) = c("trial","population")

  sd.tab = expanded.data %>%
    dplyr::group_by(trial == 0) %>%
    dplyr::summarise_all(sd) %>%
    select(-trial,-`trial == 0`,-`X.Intercept.`) %>%
    t() %>%
    as.data.frame()

  names(sd.tab) = c("trial_sd","population_sd")

  covariate_table = means.tab %>%
    bind_cols(sd.tab) %>%
    mutate(ASMD = round(abs((trial - population)/population_sd),3)) %>%
    select(trial, population, ASMD)

  row.names(covariate_table) = setdiff(names(expanded.data),c("trial","X.Intercept."))

  return(covariate_table)
}
