#'@title FitMID.
#'@description \code{FitMID} will ...
#'@param md Measured distribution. Normalized (i.e. sum=1) raw intensity vector.
#'@param td Theoretical intensity distribution (using function 'CalcTheoreticalMDV').
#'@param r matrix of considere fragments and their potential occurence.
#'@param mid_fix May provide a numeric vector used as a given MID. Allows to estimate \code{r} individually.
#'@param prec Precision of the estimation of MID, set to 1\% as default.
#'@param trace_steps For testing purposes. Print the results of intermediate steps to console.
#'@param penalize Numeric exponent penalizing solutions with low M+H occurence. Formula is 1+3*(1-x)^penalty. NA to omit penalizing.
#'@importFrom stats median
#'@importFrom utils head
#'@importFrom plyr ldply alply
#'@return Fitted MID with attributes.
#'@examples
#'\dontrun{
#'int <- c(1560, 119203, 41927, 16932, 4438)
#'names(int) <- c(-2, 0, 1, 2, 3)
#'fml <- "C19H37NO4Si3"
#'td <- CalcTheoreticalMDV(fml=fml)
#'known_frags <- unlist(list("M+H"=0,"M+"=-1,"M-H"=-2,"M+H2O-CH4"=+2))
#'r <- matrix(rep(c(0,1), 4), nrow=2, dimnames=list(NULL, names(known_frags)))
#'FitMID(md=int, td=td, trace=TRUE, frag=c(-2, 0, 1, 2, 3), r=r)
#'}
#'@keywords internal
FitMID <- function(md=NULL, td=NULL, r=NULL, mid_fix=NULL, prec=0.01, trace_steps=FALSE, penalize=NA) {

  # potential parameters
  step_increment = 0.5 # by how much is step reduced every time d2=d1*step_increment
  known_frags <- unlist(list("M+H"=0,"M+"=-1,"M-H"=-2,"M+H2O-CH4"=+2))
  if (all(r[1,]==0) & all(r[2,]==1)) limits <- NULL else limits <- r

  # set up r_fixed for internal use
  r_fixed <- ifelse(all(apply(r,2,diff)==0), TRUE, FALSE)

  # default return value
  if (sum(md)==0) {
    out <- rep(NA, nrow(td))
    names(out) <- row.names(td)
    attr(out, "err") <- unlist(list("err"=NA))
    if (prod(dim(r))>1) attr(out, "ratio") <- apply(r,2,stats::median)/sum(apply(r,2,stats::median)) else attr(out, "ratio") <- r
    attr(out, "ratio_status") <- ifelse(r_fixed, "fixed", "estimated")
    attr(out, "mid_status") <- ifelse(!is.null(mid_fix), "fixed", "estimated")
    return(out)
  }

# establish starting parameters and stepsize...
  # ...for mid
  if (!is.null(mid_fix)) {
    mid_start <- mid_fix; names(mid_start) <- rownames(td)
    mid_steps <- 0
  } else {
    mid_start <- rep(0.5,nrow(td)); names(mid_start) <- rownames(td)
    mid_steps <- 0.5
  }
  while (min(mid_steps)>prec) mid_steps <- c(mid_steps, mid_steps[length(mid_steps)]*step_increment)

  # optimize isotopomer distribution
  for (dst in mid_steps) {
    mid_local <- poss_local(vec=mid_start, d=dst, prec = prec/10, limits=NULL, length.out=3)
    if (trace_steps) {
      cat(paste("\nTesting", nrow(mid_local), "MID solutions."))
      cat(paste("\nUsing stepwidth for MID:", round(dst,5)))
    }
    test_mid <- plyr::alply(.data=as.data.frame(mid_local), 1, function (x) {
      reconstructed_mid <- apply(td*unlist(x),2,sum)
      #if ()
      #browser()
      #ReconstructMID(mid=unlist(x), r=unlist(list("M+H"=1)), fml=fml)
      ## !!! is this reconstruction above correct ??? Dont we need to multiply by rows
      if (r_fixed) {
        r_steps <- 0.5
      } else {
        r_steps <- c(0.25,0.1,0.05,0.02,0.01)
      }
      r_start <- apply(r,2,stats::median)
      d <- 0.5
      # optimize fragment ratios
      for (rst in r_steps) {
        r_local <- poss_local(vec=r_start, d=d, prec = prec/10, limits = limits, by=rst)
        if (nrow(r_local)>=1) {
          # update d for next iteration
          d <- rst
          # test fragment ratio distributions
          test_r <- apply(r_local, 1, function (y) {
            calc_mid_error(md=md, reconstructed_mid=reconstructed_mid, best_r=y)
          })
          w_r_errs <- weight_errors(rMpH = r_local[,"M+H"], errs = test_r, penalize = penalize)
          #r_start <- r_local[which.min(w_r_errs),]
          r_start <- r_local[which.min(w_r_errs),,drop=F]
          if (nrow(r_start)==0) { warning("nrow(r_start) was empty") }
          r_start <- r_start[1,]
          if (is.null(names(r_start))) names(r_start) <- colnames(r_local)
        }
      }
      best_r <- r_start
      mid_err <- calc_mid_error(md=md, reconstructed_mid=reconstructed_mid, best_r=best_r)
      return(list("err"=mid_err, "r"=round(best_r,4)))
    })
    #browser()
    #reconstructed_mid <- apply(td*unlist(attr(test_mid,"split_labels")["3",]),2,sum)
    #reconstructed_mid <- apply(td*unlist(attr(test_mid,"split_labels")["83",]),2,sum)
    #best_r <- test_mid[["1"]]$r
    #calc_mid_error(md=md, reconstructed_mid=reconstructed_mid, best_r=best_r)
    w_m_errs <- weight_errors(rMpH = sapply(test_mid,function(x){x$r["M+H"]}), errs = sapply(test_mid,function(x){x$err}), penalize = penalize)
    best <- which.min(w_m_errs)
    if (trace_steps) {
      # cat("\n\nObserved fragment ratio ranges:\n")
      # print(apply(sapply(test_mid, function(x){ x$r }),1,range))
      cat("\nTop candidates found:\n")
      tmp_print <- plyr::ldply(test_mid, function(x) { c("_____"="     ", formatC(x = x[["r"]], format = "f", digits = 2), "_____"="     ", "error"=formatC(x[["err"]], format="e", digits=2)) })
      tmp_print[,1:length(mid_start)] <- round(100*tmp_print[,1:length(mid_start)],2)
      tmp_print <- cbind(tmp_print, w_m_errs)
      print(utils::head(tmp_print[order(w_m_errs),]))
      selected_value <- readline(prompt="Type [row_number+enter] to continue stepwise or [enter] without any number to continue to end:")
      if (selected_value=="") {
        trace_steps <- FALSE
        mid_start <- mid_local[which.min(w_m_errs),]
      } else {
        mid_start <- mid_local[as.numeric(selected_value),]
      }
    } else {
      mid_start <- mid_local[which.min(w_m_errs),]
    }
    r_fin <- test_mid[[which.min(w_m_errs)]]$r
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
