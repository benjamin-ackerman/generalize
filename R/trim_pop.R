### This function works when you have a stacked data set, trial and population, with an indicator variable "trial"
## formula: formula specifying trial membership variable and covariates used to adjust
## data: if specified, can just list S and X variable names?
## just_population: logical - TRUE: spit out just the population data, FALSE: spit out full data (false = default)

trim_pop <- function(formula, data){
  trial_membership = all.vars(formula[[2]])
  covariates = all.vars(formula[[3]])

  if(anyNA(match(names(table(data[,trial_membership])),c("0","1")))){
    stop("Sample Membership variable must be coded as `0` (not in trial) or `1` (in trial)",call. = FALSE)
  }

  trial_dat = data[which(data[,trial_membership]==1),covariates]
  if(length(covariates)==1){
    trial_dat = data.frame(trial_dat)
    names(trial_dat) = covariates
    }

  covariate_bounds = function(covariate){
    if(is.factor(trial_dat[,covariate])){

      trial_levels = levels(droplevels(trial_dat)[,covariate])
      return(rownames(trim_dat)[which(!trim_dat[,covariate] %in% trial_levels)])
    }

    if(is.numeric(trial_dat[,covariate])){
      trial_bounds = c(min(trial_dat[,covariate],na.rm=TRUE),max(trial_dat[,covariate],na.rm=TRUE))
      return(rownames(data)[which(!(data[,covariate] >= trial_bounds[1] & data[,covariate] <= trial_bounds[2]))])
    }
  }

  bound_violations = purrr::map(covariates,covariate_bounds)

  missing_rows=unique(unlist(bound_violations))

  return(data[which(!rownames(data) %in% missing_rows),])
}
