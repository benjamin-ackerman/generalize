### This function works when you have a stacked data set, trial and population, with an indicator variable "trial"

## S: binary sample membership variable (1 = trial, 0 = population/not in trial)
## X: vector, matrix or dataframe containing baseline covariates
## data: if specified, can just list S and X variable names?
## just_population: logical - TRUE: spit out just the population data, FALSE: spit out full data (false = default)

trim_pop <- function(data){
    trim_dat = data
    trim_dat$row.names = row.names(trim_dat)

    trial_dat = trim_dat[which(trim_dat$trial == 1),covariates]

    covariate_bounds = function(covariate){
      if(is.factor(trial_dat[,covariate])){
        return(levels(droplevels(trial_dat)[,covariate]))
      }
      if(is.numeric(trial_dat[,covariate])){
        return(c(min(trial_dat[,covariate],na.rm=TRUE),max(trial_dat[,covariate],na.rm=TRUE)))
      }
    }
    trial_levels = map(covariates,covariate_bounds)
    names(trial_levels) = covariates

    missing_factors=list()
    for(i in names(trial_dat)){
      missing_factors[[which(names(trial_dat)==i)]] = trim_dat[which(!trim_dat[,i] %in% trial_levels[[i]]),"row.names"]
    }

    missing_rows=unique(unlist(missing_factors))
    trim_dat[which(!trim_dat$row.names %in% missing_rows),]
}
