#'@title ReconstructMID.
#'@description \code{ReconstructMID} will reconstruct a measured GC-APCI-MS spectrum
#'  of a compound given its true MID and the fragment ratio.
#'@details tbd
#'@param mid A numeric vector with sum=1 and length of C atoms +1.
#'@param r Fragment ratios. A numeric vector with sum=1.
#'@param fml A compound formula.
#'@param cutoff Remove values below this threshold from output vector.
#'@return A reconstructed MID.
#'@keywords internal
#'@examples
#'fml <- "C9H20O3Si2"
#'attr(fml,"nmz") <- 9
#'attr(fml,"nbio") <- 3
#'# make up some fake measurement data for this based on assumed 10% labelling at M3
#'mid <- c(0.9,0,0,0.1)
#'r <- unlist(list("M+H"=0.8, "M+"=0.1, "M+H2O-CH4"=0.1))
#'(rMID <- CorMID:::ReconstructMID(mid=mid, r=r, fml=fml, cutoff=0.001))
#'CorMID:::plot_rawMID(rMID)
#'\dontrun{
#'CorMID::CorMID(int = rMID, fml=fml, prec=0.001, r=r, trace_steps = TRUE)
#'}
ReconstructMID <- function(mid=NULL, r=unlist(list("M+H"=1)), fml=NULL, cutoff=0.01) {
  # ensure that sum(mid)==1
  mid <- mid/sum(mid)
  if (is.null(attr(fml, "nbio"))) attr(fml, "nbio") <- length(mid)-1
  if (is.null(attr(fml, "nmz"))) attr(fml, "nmz") <- length(mid)+1
  td <- CalcTheoreticalMDV(fml=fml)
  n <- ncol(td)
  s <- colSums(td*mid)
  # ensure that only (and all) specified r are available
  frag <- unlist(list("M-H"=0,"M+"=0,"M+H"=0,"M+H2O-CH4"=0))
  r <- r[names(r) %in% names(frag)]
  frag[names(r)] <- r
  r <- frag; r <- r/sum(r)
  # reconstruct rawMID
  s <- r["M-H"]*s + c(0, r["M+"]*s[1:(n-1)]) + c(0, 0, r["M+H"]*s[1:(n-2)]) + c(rep(0,4), r["M+H2O-CH4"]*s[1:(n-4)])
  names(s) <- paste0("M", formatC(-2:(length(s)-3), flag="+"))
  if (is.finite(cutoff)) s <- s[s>=cutoff]
  return(s)
}
