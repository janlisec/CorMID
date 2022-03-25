#'@title weight_errors.
#'@description \code{weight_errors} will apply a weighting function to a numeric vector of errors.
#'@param rMpH A numeric vector of fragment ratios for the M+H fragment.
#'@param errs A numeric vector of errors.
#'@param penalize A single numeric or NA to omit weighting.
#'@param max_panel Constant to modify weighting function
#'@return A numeric vector of length(x).
#'@keywords internal
#'@examples
#'rMpH <- seq(0,1,0.1)
#'errs <- rep(1,11)
#'CorMID:::weight_errors(rMpH=rMpH, errs=errs, penalize=6)
#'plot(1,1,xlim=c(0,1), ylim=c(0,11), type="n", xlab="rMpH", ylab="Weighted Error")
#'for (i in 1:10) {
#'lines(x=rMpH, y=CorMID:::weight_errors(rMpH=rMpH, errs=errs, penalize=i), col=i)
#'}
#'legend(x="topright", legend=1:10, fill=1:10, title="penalize")
weight_errors <- function(rMpH=NULL, errs=NULL, penalize=NA, max_panel=9) {
  if (is.finite(penalize)) {
    fac <- 1+max_panel*(1-rMpH)^penalize
    return(errs*fac)
  } else {
    return(errs)
  }
}
