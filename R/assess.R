#' Assess Generalizability of Randomized Trial
#'
#' @param formula an object of class "formula". The formula specifying the model for trial participation.  Lefthand side should be a binary variable indicating trial membership, and righthand side should contain pre-treatment covariates measured in data set.
#' @param data a data frame containing the variables specified in the model
#' @param method choose method to predict the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are random forests ("rf") and lasso ("lasso")
#' @param trim.pop logical, if TRUE - trim target population so that covariates do not exceed bounds of trial covariates
#' @param is.data.disjoint logical, if TRUE - trial and population data are considered independent.  If FALSE - trial is proper subset of population data.  See details for implications on weighting methods.
#' @return
#' @examples
#' gen_weights(trial ~ age + sex + race, data = ctn_data)
#' gen_weights(trial ~ age + sex + race, data = ctn_data, method = 'rf')

assess = function(formula, data, method = "lr", trim.pop = TRUE, is.data.disjoint = TRUE){
  ### Get variable names from the formula ###
  trial_membership = all.vars(formula)[1]
  covariates = all.vars(formula)[-1]

  ### Checks ###

  # if (missing(data)) {
  #   data <- environment(formula)
  #   #stop("Data must be specified.", call. = FALSE)
  # }

  if(class(formula) != "formula"){
    stop("Must enter a valid formula!",call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if(anyNA(match(all.vars(formula),names(data)))){
    missing_variables = all.vars(formula)[is.na(match(all.vars(formula),names(data)))]
    stop(paste0(paste(missing_variables,collapse = ", ")," not in the data provided"))
  }

  if(!length(unique(data[,trial_membership])) == 2){
    stop("Trial Membership variable not binary", call. = FALSE)
  }

  ### Clean up data from missing values ###
  data = data[rownames(na.omit(data[,all.vars(formula)])),]

  if(trim.pop = TRUE){
    trimmed_data = trim_pop(formula, data = data)$trimmed_data

    ### Find number of population members excluded from trimming
    n_excluded = trim_pop(formula, data = data)$n_excluded

    ### Calculate participation probabilities, weights, generalizability index
    gen_weights_object = gen_weights(formula, data = trimmed_data, method = method, is.data.disjoint = is.data.disjoint)
    participation_probs = gen_weights_object$participation_probs
    weights = gen_weights_object$weights

    g_index = gen_index(participation_probs$probs_trial, participation_probs$probs_population)

    return(list(
      n_excluded = n_excluded,
      g_index = g_index
    ))
  }

  if(trim.pop = FALSE){
    ### Calculate participation probabilities, weights, generalizability index
    gen_weights_object = gen_weights(formula, data = data, method = method, is.data.disjoint = is.data.disjoint)
    participation_probs = gen_weights_object$participation_probs
    weights = gen_weights_object$weights

    g_index = gen_index(participation_probs$probs_trial, participation_probs$probs_population)

    return(list(
      g_index = g_index
    ))
  }
}
