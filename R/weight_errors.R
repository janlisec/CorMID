#' @title weight_errors.
#' @description \code{weight_errors} will apply a weighting function to a numeric vector of errors.
#' @param rMpH A numeric vector of fragment ratios for the M+H fragment.
#' @param errs A numeric vector of errors.
#' @param penalize A single numeric or NA to omit weighting.
#' @param max_panel Constant to modify weighting function
#' @return A numeric vector of length(x).
#' @keywords internal
#' @noRd
weight_errors <- function(rMpH=NULL, errs=NULL, penalize=NA, max_panel=9) {
  if (is.finite(penalize)) {
    fac <- 1+max_panel*(1-rMpH)^penalize
    return(errs*fac)
  } else {
    return(errs)
  }
}
