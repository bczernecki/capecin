#' Export all missing functions
#'
#' TODO
#' To be completed..
#' @export
#' @noRd

getthe = function(p, t, td, q) {

  if ((td - t) > -0.1) {
    tlcl = t
  } else {
  tlcl = 56.0 + ((td - 56.0)^(-1) + 0.00125*log(t/td) )^(-1)
  }

 getthe = t*((100000.0/p)^(0.2854*(1.0 - 0.28*q)))*exp(((3376.0/tlcl) - 2.54)*q*(1.0 + 0.81*q))
 return(getthe)

}

#getthe(p = 1000, t = 15, td = 10, q = 0.001)


getqvs = function(p, t) {
  eps = 287.04/461.5
  es = 611.2*exp(17.67*(t - 273.15)/(t - 29.65))
  getqvs = eps*es/(p - es)
  return(getqvs)
}
#getqvs(1000, 15)

getqvi = function(p, t) {
  eps = 287.04/461.5
  es = 611.2*exp(21.8745584*(t - 273.15)/(t - 7.66))
  getqvi = eps*es/(p - es)
  return(getqvi)
}

