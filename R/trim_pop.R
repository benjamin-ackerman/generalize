#' Trim target population covariates to be within bounds of trial
#'
#' @param formula an object of class "formula". The formula specifying the model for trial participation.  Lefthand side should be a binary variable indicating trial membership, and righthand side should contain pre-treatment covariates measured in data set.
#' @param data a data frame containing the variables specified in the model
#' @param just_population logical. if FALSE (default), function returns full data set, if TRUE, function only returns trimmed population data.
#' @return \code{trim_pop} returns a data frame, where the target population covariates do not exceed the bounds of the trial covariates
#' @examples
#' trim_pop(trial ~ age + sex + race, data = ctn_data)

trim_pop <- function(formula, data, just_population = FALSE){

  trial_membership = all.vars(formula[[2]])
  covariates = all.vars(formula[[3]])

  if(class(formula) != "formula"){
    stop("Must enter a valid formula!",call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if(anyNA(match(all.vars(formula),names(data)))){
    missing_variables = all.vars(formula)[is.na(match(all.vars(formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," are not variables in the data provided"))
  }

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
      return(rownames(data)[which(!data[,covariate] %in% trial_levels)])
    }

    if(is.numeric(trial_dat[,covariate])){
      trial_bounds = c(min(trial_dat[,covariate],na.rm=TRUE),max(trial_dat[,covariate],na.rm=TRUE))
      return(rownames(data)[which(!(data[,covariate] >= trial_bounds[1] & data[,covariate] <= trial_bounds[2]))])
    }
  }

  bound_violations = purrr::map(covariates,covariate_bounds)

  missing_rows=unique(unlist(bound_violations))

  trimmed_data = data[which(!rownames(data) %in% missing_rows),]

  if(just_population == TRUE){
    return(trimmed_data[which(trimmed_data[,trial_membership] == 0),])
  }

  if(just_population == FALSE){
    return(trimmed_data)
  }

}
