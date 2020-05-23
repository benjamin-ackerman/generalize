#' Create Covariate Balance Table
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param weighted_table should the output be a weighted table?
#' If \code{TRUE}, then \code{\link{weighting}}
#' @param selection_method method to estimate the probability of trial
#' participation.  Default is logistic regression ("lr").
#' Other methods supported are Random Forests ("rf") and Lasso ("lasso"),
#' passed to \code{\link{weighting}}
#' @param sl_library vector of SuperLearner library methods. If `selection_method` = 'super', specify names of methods to include in library. Default is NULL.
#' @param survey_weights variable name of population data's complex survey weights. Default is \code{FALSE}: if \code{FALSE}, then population data do not come a complex survey and weights do not need to be incorporated in estimation.
#' @param trim_weights logical. If \code{TRUE}, then trim the weights to the value specified in `trim_pctile`. Default is \code{FALSE}.
#' @param trim_pctile numeric. If `trim_weights` is \code{TRUE}, then specify what percentile weights should be trimmed to. Default is 0.97.
#' @param is_data_disjoint logical. If \code{TRUE}, then trial and population data
#'  are considered independent.  This affects calculation of the weights -
#'  see details for more information.
#' @importFrom dplyr funs

#' @export
covariate_table <- function(trial, selection_covariates, data,
                            weighted_table = FALSE,
                            selection_method = "lr",
                            sl_library = NULL,
                            survey_weights = FALSE,
                            trim_weights=FALSE,
                            trim_pctile = .97,
                            is_data_disjoint = TRUE){
  . = population = pooled_sd = ASMD = NULL
  rm(list = c( "population", "pooled_sd", "ASMD", "."))

  data = data %>%
    tidyr::drop_na(selection_covariates) %>%
    as.data.frame()

  if(survey_weights == FALSE){
    data$s_weights = 1
  } else{
    data$s_weights = ifelse(data[,trial] == 1, 1, data[,survey_weights])

    ### NORMALIZE THE SURVEY WEIGHTS TO SUM TO THE NUMBER OF SURVEY PARTICIPANTS:
    normalize_factor = mean(data$s_weights[which(data[,trial] == 0)], na.rm = TRUE)
    data$s_weights[which(data[,trial] == 0)] = data$s_weights[which(data[,trial] == 0)]/normalize_factor
  }

  if(weighted_table == TRUE){
    data$weights = weighting(outcome = NULL, treatment = NULL, trial = trial,
                             selection_covariates = selection_covariates, data = data,
                             selection_method = selection_method,
                             sl_library = sl_library,
                             survey_weights=survey_weights,
                             trim_weights=trim_weights,
                             trim_pctile=trim_pctile,
                             is_data_disjoint = is_data_disjoint)$weights

    data$weights = ifelse(data[,trial] == 0, data$s_weights, data$weights)
  } else {
    data$weights = data$s_weights
  }

  dmy <- caret::dummyVars(" ~ .", data = data[,selection_covariates], sep="_")
  expanded.data = data.frame(trial = data[, c(trial)],
                             predict(dmy, newdata = data[,selection_covariates]),
                             weights = data$weights)
  names(expanded.data)=stringr::str_replace_all(names(expanded.data),"\\.","-")

  means.tab = expanded.data %>%
    dplyr::group_by(trial) %>%
    dplyr::summarise_at(names(expanded.data)[-1], dplyr::funs(weighted.mean(., weights))) %>%
    dplyr::select(-`weights`) %>%
    t() %>%
    as.data.frame()

  names(means.tab)[which(means.tab["trial",]=="1")] = "trial"
  names(means.tab)[which(means.tab["trial",]=="0")] = "population"

  means.tab = means.tab[!rownames(means.tab) %in% c("trial", "X.Intercept."),]

  n_trial = as.numeric(table(expanded.data$trial)["1"])
  n_pop = as.numeric(table(expanded.data$trial)["0"])

  # as per
  # https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weighvar.pdf
  sd.tab = expanded.data %>%
    dplyr::group_by(trial) %>%
    dplyr::summarise_at(names(expanded.data)[which(!names(expanded.data)%in% c("trial","X.Intercept."))],
                        dplyr::funs(sum(weights * (. - weighted.mean(.,weights))^2)/
                                      ((length(weights) - 1) * sum(weights) / length(weights))
                                    )) %>%
    dplyr::select(-`weights`) %>%
    t() %>%
    as.data.frame()

  sd.tab$pooled_sd = sqrt(((n_trial - 1)*sd.tab$V1 + (n_pop - 1)*sd.tab$V2)/(n_trial + n_pop - 2))

  names(sd.tab)[which(sd.tab["trial",]=="1")] = "trial_var"
  names(sd.tab)[which(sd.tab["trial",]=="0")] = "population_var"
  names(sd.tab)[3] = "pooled_sd"

  sd.tab = sd.tab[!rownames(sd.tab) %in% c("trial", "X.Intercept."),]

  covariate_table = means.tab %>%
    dplyr::bind_cols(sd.tab) %>%
    dplyr::mutate(ASMD = round(abs((trial - population)/pooled_sd),3)) %>%
    dplyr::select(trial, population, ASMD)

  if(survey_weights != FALSE){
    names(covariate_table)[which(names(covariate_table) == "population")] = "population (weighted)"
  }

  if(weighted_table == FALSE){
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial","X.Intercept.","weights"))
  }

  if(weighted_table == TRUE){
    names(covariate_table)[which(names(covariate_table) == "trial")] = "trial (weighted)"
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial","X.Intercept.","weights"))
  }

  return(covariate_table)
}
