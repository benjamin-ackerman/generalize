gen_index <- function(dat1B,dat2B) {
  # bandwidth
  h = function(x){
    n = length(x)
    return((4*sqrt(var(x))^5/(3*n))^(1/5))
  }#kernel density
  kg = function(x,data){
    hb = h(data) #bin width
    k = r = length(x)
    for(i in 1:k)
      r[i] = mean(dnorm((x[i]-data)/hb))/hb
    return(r)
  }##B index calculation
  return( as.numeric(integrate(function(x) sqrt(kg(x,dat1B)*kg(x,dat2B)),-Inf,Inf)$value))
}
