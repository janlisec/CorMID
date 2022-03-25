#'@title calc_mid_error.
#'@description \code{calc_mid_error} will compute the error of a theoretical and a estimated mid.
#'@param md Normalized measured intensities
#'@param reconstructed_mid A reconstructed MID based on a true MID and a theoretical distribution.
#'@param best_r A named numeric vector of fragment ratios.
#'@return A numeric vector of length(x).
#'@keywords internal
#'@examples
#'\dontrun{
#'md <- c(50,100,10,0,50,50); md <- md/sum(md); names(md) <- -1:4
#'mid <- c(0.7,0,0,0,0.3); names(mid) <- paste0("M",0:4)
#'reconstructed_mid <- apply(t(CalcTheoreticalMDV(fml="C19H37NO4Si3")[1:5,1:6])*mid,1,sum)
#'best_r <- unlist(list("M+H"=0.9,"M+"=0.1))
#'CorMID:::calc_mid_error(md=md, reconstructed_mid=reconstructed_mid, best_r=best_r)
#'}
calc_mid_error <- function(md=NULL, reconstructed_mid=NULL, best_r=NULL) {
  known_frags <- unlist(list("M+H"=0,"M+"=-1,"M-H"=-2,"M+H2O-CH4"=+2))
  frag <- as.numeric(gsub("M","",names(md)))
  out <- rep(0, length(reconstructed_mid))
  names(out) <- names(reconstructed_mid)
  for (i in 1:length(best_r)) {
    # compute necessary trailing zeros
    n0 <- rep(0, which(frag==known_frags[names(known_frags) %in% names(best_r)[i]])-1)
    # reconstruct the MID for this adduct
    rMID <- reconstructed_mid*unlist(best_r[i])
    # add the MID for this adduct
    out <- out + c(n0, rMID)[1:length(out)]
  }
  # normalized and calculate error for this mid
  # $$ the following line led to serious problems and was outcommented on 2022-03-16 by JL
  #if (sum(out)!=1 & sum(out)!=0) out <- out/sum(out)
  mid_err <- sqrt(sum((out-md)^2))
  return(mid_err)
}
