#' Methamphetamine Trial and Target Population Data
#'
#' A data set containing data from CSP-1025, a randomized trial from the CTN database, and TEDS-A 2014, a nationally-representative collection of admissions records to publicly-funded substance use facilities.
#'
#' @format A data frame with 135404 rows and 15 variables:
#' \describe{
#'   \item{ID}{subject ID}
#'   \item{trial}{trial membership indicator: 1 = belongs to CSP-1025, 0 = belongs to TEDS-A-2014}
#'   \item{age}{age category: 12-14, 15-17, 18-20, 21-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55+}
#'   \item{sex}{sex: male, female}
#'   \item{race}{race: Asian, Black, Native Hawaiian, White, Other}
#'   \item{ethnicity}{ethnicity: Hispanic/Latino, Not Hispanic/Latino, Unknown/Not Given}
#'   \item{marital_status}{marital status: Divorced/Not Married, Married/Partnered}
#'   \item{education}{years of education: <9, 9-11, 12, 13-15, 16+}
#'   \item{employment}{employment status: Full time, Part time, Not in labor force, Unemployed}
#'   \item{methprior}{Was methamphetamine used in the prior 7 days? 1 = yes, 0 = no}
#'   \item{treat}{(only measured in trial) treatment indicator: 1 = treatment, 0 = placebo}
#'   \item{studyretention}{(only measured in trial) Did the participant remain in the trial until the end? 1 = yes, 0 = no}
#'   \item{methfollowup}{(only measured in trial) Did the participant report using methamphetamine in trial followup? 1 = yes, 0 = no}
#'   \item{methwk11}{(only measured in trial) Did the participant have 3 methamphetamine-free drug screenings in week 11 of the trial? 1 = yes, 0 = no}
#'   \item{methwk12}{(only measured in trial) Did the participant have 3 methamphetamine-free drug screenings in week 12 of the trial? 1 = yes, 0 = no}
#' }
#' @source CTN data: \url{https://clinicaltrials.gov/ct2/show/study/NCT00345371}
#' @source TEDS-A data: \url{https://www.datafiles.samhsa.gov/study-dataset/treatment-episode-data-set-admissions-2014-teds-2014-ds0001-nid16950}

"meth_data"
