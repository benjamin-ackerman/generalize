##### CITE BETH TIPTON'S PAPER IN THE DOCUMENTATION #####

#' Calculate `generalizability index` to describe how similar or different trial is from population
#'
#' @param dat1B vector of probabilities of trial participation among individuals in the trial
#' @param dat2B vector of probabilities of trial participation among individuals in the population
#' @return the generalizability index, a value between 0 and 1, where scores greater than 1 indicate greater similarity (see Tipton paper for description)

#' @export
gen_index <- function(dat1B,dat2B) {
  #kernel density
  kg = function(x,data){
    # Bandwidth
    n = length(data)
    hb = (4*sqrt(var(data))^5/(3*n))^(1/5)

    # Kernel density
    k = r = length(x)
    for(i in 1:k)
      r[i] = mean(dnorm((x[i]-data)/hb))/hb
    return(r)
  }

  ##B index calculation
  return( as.numeric(integrate(function(x) sqrt(kg(x,dat1B)*kg(x,dat2B)),-Inf,Inf)$value))
}
