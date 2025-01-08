#' @title FitMID.
#' @description \code{FitMID} will ...
#' @param md Measured distribution. Normalized (i.e. sum=1) raw intensity vector.
#' @param td Theoretical intensity distribution (using function 'CalcTheoreticalMDV').
#' @param r matrix of consider fragments and their potential occurrence.
#' @param mid_fix May provide a numeric vector used as a given MID. Allows to estimate \code{r} individually.
#' @param prec Precision of the estimation of MID, set to 1\% as default.
#' @param trace_steps For testing purposes. Print the results of intermediate steps to console.
#' @param penalize Numeric exponent penalizing solutions with low M+H occurrence. Formula is 1+3*(1-x)^penalty. NA to omit penalizing.
#' @return Fitted MID with attributes.
#' @keywords internal
#' @noRd
#' @importFrom Rcpp sourceCpp
#' @useDynLib CorMID, .registration = TRUE
FitMID <- function(md=NULL, td=NULL, r=NULL, mid_fix=NULL, prec=0.01, trace_steps=FALSE, penalize=NA) {

  # potential parameters
  step_increment = 0.5 # by how much is step reduced every time d2=d1*step_increment
  known_frags <- CorMID::known_frags
  if (prod(dim(r))>1 && all(r[1,]==0) & all(r[2,]==1)) limits <- NULL else limits <- r

  # default return value
  if (sum(md)==0) {
    message("No finite intensity values provided. Return NA vector.")
    out <- as.numeric(rep(NA, ifelse(is.null(td), length(md), nrow(td))))
    if (is.null(td)) names(out) <- names(md) else names(out) <- row.names(td)
    attr(out, "err") <- unlist(list("err"=NA))
    if (prod(dim(r))>1) attr(out, "ratio") <- apply(r,2,stats::median)/sum(apply(r,2,stats::median)) else attr(out, "ratio") <- r
    attr(out, "ratio_status") <- ifelse(prod(dim(r))>1 && all(apply(r,2,diff)==0), "fixed", "estimated")
    attr(out, "mid_status") <- ifelse(!is.null(mid_fix), "fixed", "estimated")
    return(out)
  } else {
    # ensure that intensity vector is normalized to sum
    md <- md/sum(md)
  }

  # set up r_fixed for internal use
  r_fixed <- ifelse(all(apply(r, 2, diff)==0), TRUE, FALSE)

  # establish starting parameters and step size...
  # ...for mid
  if (!is.null(mid_fix)) {
    mid_start <- stats::setNames(mid_fix, rownames(td))
    mid_steps <- 0
  } else {
    mid_start <- stats::setNames(rep(0.5,nrow(td)), rownames(td))
    mid_steps <- 0.5
    #mid_steps <- 0.5*step_increment^(0:10)
  }
  while (min(mid_steps)>prec) mid_steps <- c(mid_steps, mid_steps[length(mid_steps)]*step_increment)

  # optimize isotopologue distribution
  # JL$$ prepare as much data as possible for error calculation function
  frag <- as.numeric(substr(names(md),2,4))
  n_md <- length(md)
  L0 <- sapply(colnames(r), function(x) { rep(0, abs(min(frag)-known_frags[x])) }, simplify = FALSE)
  shift_r <- sapply(colnames(r), function(x) { abs(min(frag)-known_frags[x]) }, USE.NAMES = FALSE)
  # $$JL

  for (dst in mid_steps) {
    #mid_local <- poss_local(vec=mid_start, d=dst, prec = prec/10, limits=NULL, length.out=3)
    mid_local <- poss_local_C(vec=mid_start, d=dst, prec = prec, length=3)
    if (trace_steps) {
      cat(paste("\nTesting", nrow(mid_local), "MID solutions."))
      cat(paste("\nUsing stepwidth [%] for MID:", round(100*dst,3)))
    }
    #browser()
    test_mid <- apply(as.data.frame(mid_local), 1, function (x) {
      pre_mid <- colSums(td*unlist(x))
      if (r_fixed) {
        r_steps <- 0.5
      } else {
        # JL 20241120 r_steps has been replaced
        #r_steps <- c(0.25,0.1,0.05,0.02,0.01)
        r_steps <- 0.25/2^(0:5)
      }
      r_start <- apply(r, 2, stats::median)
      d <- 0.5
      # optimize fragment ratios
      for (rst in r_steps) {
        #r_local <- poss_local(vec=r_start, d=d, prec = prec/10, limits = limits, by=rst)
        #browser()
        r_local <- poss_local_C(vec=r_start, d=d, prec = prec/10, limits = limits, by=rst)
        if (FALSE) {
          r_local_old <- poss_local(vec=r_start, d=d, prec = prec/10, limits = limits, by=rst)
          if (!identical(r_local, r_local_old)) browser()
        }
        if (nrow(r_local)>=1) {
          # update d for next iteration
          d <- rst
          # test fragment ratio distributions
          test_r <- apply(r_local, 1, function (y) {
            #calc_mid_error(md=md, reconstructed_mid=pre_mid, best_r=y, n_md=n_md, L0=L0)
            calc_mid_error_C(md, pre_mid, y, shift_r)
          })
          w_r_errs <- weight_errors(rMpH = r_local[,"M+H"], errs = test_r, penalize = penalize)
          r_start <- r_local[which.min(w_r_errs)[1],,drop=F]
          if (nrow(r_start)==0) { warning("nrow(r_start) was empty") }
          r_start <- r_start[1,]
          if (is.null(names(r_start))) names(r_start) <- colnames(r_local)
        }
      }
      best_r <- r_start
      #mid_err <- calc_mid_error(md=md, reconstructed_mid=pre_mid, best_r=best_r, n_md=n_md, L0=L0)
      mid_err <- calc_mid_error_C(md, pre_mid, best_r, shift_r)
      return(list("err"=mid_err, "r"=round(best_r,4)))
    })
    # compute weighted errors (penalty for solutions with low M+H ion, which is unlikely)
    w_m_errs <- weight_errors(
      rMpH = sapply(test_mid, function(x) { x$r["M+H"] }),
      errs = sapply(test_mid, function(x) { x$err }),
      penalize = penalize
    )
    # get new starting position of MID either through user input or as optimal solution of previous step via min of weighted errors
    if (trace_steps) {
      #browser()
      new_mid_start <- select_mid_start_manual(test_mid, mid_local, w_m_errs)
      if (is.null(new_mid_start)) {
        trace_steps <- FALSE
        new_mid_start <- names(which.min(w_m_errs))
      }
    } else {
      new_mid_start <- names(which.min(w_m_errs))
    }
    mid_start <- mid_local[new_mid_start,]
    r_fin <- test_mid[[new_mid_start]]$r
  }

  # ensure once more that sum=100 and result is given in %
  mid_fin <- 100*unlist(mid_start)/sum(mid_start)
  attr(mid_fin, "err") <- min(sapply(test_mid,function(x){x$err}))
  names(attr(mid_fin, "err")) <- "err"
  attr(mid_fin, "ratio") <- r_fin
  attr(mid_fin, "ratio_status") <- ifelse(r_fixed, "fixed", "estimated")
  attr(mid_fin, "mid_status") <- ifelse(!is.null(mid_fix), "fixed", "estimated")
  return(mid_fin)
}
