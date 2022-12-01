#' @title poss_local.
#' @description \code{poss_local} will compute a matrix of possibilities.
#' @details Within the approximation process we need to check various hypotheses
#'  of MID and r combinations. A non-redundant set of posible combinations can
#'  be computed with this \code{poss_local}.
#' @param vec The starting vector sum(vec) should be 1.
#' @param d The maximum allowed deviation for each element in vec.
#' @param prec recision of allowed errors.
#' @param limits A 2-row matrix with lower and upper boundaries for the result vectors.
#' @param ... Passed to function \code{seq}. Either by or length.out (see examples in test-poss_local.R).
#' @return A matrix with rowSums\~1 and within the limits defined by vec and d.
#' @keywords internal
#' @noRd
poss_local <- function(vec=NULL, d=NULL, prec=0.001, limits=NULL, ...) {
  # helper
  fnc_lim <- function(lp=NULL, r=r) {
    flt <- sapply(1:ncol(lp), function(i) {
      lp[,i]>=(r[1,i]-prec) & lp[,i]<=(r[2,i]+prec)
    })
    if (is.null(nrow(flt))) {
      flt <- all(flt)
    } else {
      flt <- apply(flt, 1, all)
    }
    out <- lp[flt,,drop=FALSE]
    return(out)
  }

  if (d==0) {
    lst <- as.list(vec)
  } else {
    lst <- lapply(vec, function(x) {seq(from = max(0, x-d), to = min(1, x+d), ...)})
  }
  tmp <- expand.grid(lst)
  # sum is often not ident(1), hence the solution to compare with precision
  flt <- abs(rowSums(tmp)-1)<=prec
  if (any(flt)) {
    lp <- tmp[flt,,drop=F]
    if (!is.null(limits)) {
      lp <- fnc_lim(lp = lp, r = limits)
    }
  } else {
    vec <- vec/sum(vec)
    lp <- t(unlist(vec))
  }
  # if (nrow(lp)==0) {
  #   lp <- t(as.matrix(vec/sum(vec)))
  # }
  return(as.matrix(lp))
}
