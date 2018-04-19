## Output: something that documents how many excluded, summarizes the "table 1" of covariate means for included/excluded
## Continuous variables: allow for some amount of "wiggle room" - threshold is not min or max

#' Subset Population so Population Covariates are within bounds of Trial Covariates
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @return \code{trim_pop} returns a data frame, where the target population covariates do not exceed the bounds of the trial covariates

trim_pop <- function(trial, selection_covariates, data){

  ##### CHECKS #####
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  if(anyNA(match(selection_covariates,names(data)))){
    stop("Not all covariates listed are variables in the data provided!",call. = FALSE)
  }

  if(!length(na.omit(unique(data[,trial]))) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  if(anyNA(match(names(table(data[,trial])),c("0","1")))){
    stop("Sample Membership variable must be coded as `0` (not in trial) or `1` (in trial)",call. = FALSE)
  }

  ##### subset trial data covariates #####
  trial_dat = data[which(data[,trial]==1),selection_covariates]

  if(length(selection_covariates)==1){
    trial_dat = data.frame(trial_dat)
    names(trial_dat) = selection_covariates
    }

  ##### find covariate bounds in the trial #####
  covariate_bounds = function(covariate){
    if(is.factor(trial_dat[,covariate])){

      trial_levels = levels(droplevels(trial_dat)[,covariate])
      return(rownames(data)[which(!data[,covariate] %in% trial_levels)])
    }

    if(is.numeric(trial_dat[,covariate])){
      trial_bounds = c(min(trial_dat[,covariate],na.rm=TRUE),max(trial_dat[,covariate],na.rm=TRUE))
      return(rownames(data)[which(!(data[,covariate] >= trial_bounds[1] & data[,covariate] <= trial_bounds[2]))])
    }
  }

  ##### find rows of population data that violate bounds #####
  bound_violations = purrr::map(selection_covariates,covariate_bounds)

  missing_rows=unique(unlist(bound_violations))

  trimmed_data = data[which(!rownames(data) %in% missing_rows),]

  ##### get rid of unused levels from factors #####
  trimmed_data = droplevels(trimmed_data)

  ##### number of rows in population data excluded #####
  n_excluded = nrow(data[which(data[,trial]==0),]) - nrow(trimmed_data[which(trimmed_data[,trial]==0),])

  out = list(n_excluded = n_excluded,
             trimmed_data = trimmed_data,
             untrimmed_data = data)

  return(out)
}
