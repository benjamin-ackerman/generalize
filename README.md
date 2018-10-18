# generalize: an R package for estimating population effects from randomized trial data

The "generalize" package contains two core functions: `assess` and `generalize`.

## `assess`
`Assess`  evaluates similarities and differences between the trial sample and the target population based on a specified list of common covariates.  This is done in a few ways: first, `assess` provides a summary table of covariate means in the trial and the population, along with absolute standardized mean differences (ASMD) between the two sources of data.  Second, `assess` estimates the probability of trial participation based on a vector of covariate names and specified method, and summarizes their distribution across the trial and target populations.  For this, logistic regression is the default method, but estimation using Random Forests or Lasso are currently supported by the package as well.  Third, `assess` utilizes the estimated trial participation probabilities to calculate the Tipton generalizability index.  Lastly, `assess` allows researchers to "trim" the target population data set so that the covariate bounds do not exceed those of the trial covariates.  This checks for any violations of the coverage assumption (that the distribution of the covariates in the population are within the bounds of the covariate distributions in the trial), and reports how many individuals in the population would be excluded.

```
assess(trial = "trial", selection_covariates = covariates, 
       data = meth_data, selection_method = "rf", trim_pop = TRUE)
```

The summary of an object created with `assess` returns the selection model, the distribution of the trial participation probabilities by data source, and the method of trial participation probability estimation.  It also contains the calculated Tipton generalizability index, the number of individuals excluded due to coverage violations, and a table of the covariate distributions.

## `generalize`
After assessing the generalizability, the `generalize` function can be used to implement methods to estimate the target population average treatment effect (TATE).  Weighting by the odds using logistic regression is the default method, though weights based on other models (Lasso or Random Forests) or using BART or TMLE are available for use as well.

```
generalize(outcome = "methfollowup", treatment = "treat", trial = "trial", 
           selection_covariates = covariates, data = meth_data, 
           method = "weighting", selection_method = "rf", trim_pop = FALSE)
```

The summary of an object created using `generalize` returns a table with the SATE and TATE estimates, along with their standard errors and 95% confidence intervals.  When weighting is used, a covariate distribution table is printed as well, where the covariate means in the trial are weighted by the generated trial participation weights.
