#'@title recMID.
#'@description \code{recMID} will reconstruct a measured GC-APCI-MS spectrum
#'  of a compound given its true MID and the fragment ratio.
#'@details \code{recMID} is basically the inverse function to \code{CorMID}.
#'  Providing a specific chemical formula together with information regarding
#'  the true MID and r, this function will compute a vector of ion intensities
#'  which can be expected in a GC-APCI-MS analysis for this compound.
#'@param mid A numeric vector with sum=1 and length of C atoms +1.
#'@param r Fragment ratios. A numeric vector with sum=1.
#'@param fml A compound formula.
#'@param cutoff Remove values below this threshold from output vector.
#'@return A reconstructed MID.
#'@export
#'@examples
#'fml <- "C9H20O3Si2"
#'mid <- c(0.9,0,0,0.1)
#'r <- list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1)
#'(rMID <- CorMID::recMID(mid=mid, r=r, fml=fml))
#'plot(rMID)
#'plot(x = rMID, ylim=c(0,max(rMID)))
#'plot(x = rMID, xlim=c(-2,12), ylim=NULL, col=2, lwd=12, las=2, xlab="label")
#'\donttest{
#'CorMID::CorMID(int = rMID, fml=fml, prec=0.001, r=unlist(r), trace_steps = TRUE)
#'}
recMID <- function(mid=NULL, r=list("M+H"=1), fml=NULL, cutoff=0.001) {
  # ensure that sum(mid)==1
  mid <- mid/sum(mid)
  nbio <- ifelse(is.null(attr(fml, "nbio")), length(mid)-1, attr(fml, "nbio"))
  nmz <- ifelse(is.null(attr(fml, "nmz")), length(mid)+3, attr(fml, "nmz"))
  td <- CalcTheoreticalMDV(fml=fml, nbio = nbio, nmz = nmz)
  n <- ncol(td)
  s <- colSums(td*mid)
  # ensure that only (and all) specified r are available
  frag <- unlist(list("M-H"=0,"M+"=0,"M+H"=0,"M+H2O-CH4"=0))
  r <- unlist(r)
  r <- r[names(r) %in% names(frag)]
  frag[names(r)] <- r
  r <- frag; r <- r/sum(r)
  # reconstruct rawMID
  s <- r["M-H"]*s + c(0, r["M+"]*s[1:(n-1)]) + c(0, 0, r["M+H"]*s[1:(n-2)]) + c(rep(0,4), r["M+H2O-CH4"]*s[1:(n-4)])
  names(s) <- paste0("M", formatC(-2:(length(s)-3), flag="+"))
  if (is.finite(cutoff)) s <- s[s>=cutoff]
  s <- s/sum(s)
  class(s) <- "recMID"
  return(s)
}

#'@rdname recMID
#'@param x Object of class recMID.
#'@param ... Further plotting parameters.
#'@importFrom graphics axis par
#'@importFrom grDevices grey
#'@export
#'@method plot recMID
plot.recMID <- function(x, ...) {
  stopifnot(is.numeric(x))
  if (is.null(names(x))) {
    names(x) <- paste("M", 0:(length(x)-1), sep="+")
  } else {
    if (all(substr(names(x),1,1)=="M")) {

    } else {
      if (all(is.finite(as.numeric(names(x))))) {
        names(x) <- paste0("M", formatC(x = as.numeric(names(x)), format = "d", flag = "+"))
      }
    }
  }
  stopifnot(all(names(x) %in% c("M-2", "M-1", paste("M", 0:12, sep="+"))))
  argg <- c(as.list(environment()), list(...))
  if (!"xlab" %in% names(argg)) xlab <- "" else xlab <- argg$xlab
  if (!"ylab" %in% names(argg)) ylab <- "Relative intensity" else ylab <- argg$ylab
  if (!"lwd" %in% names(argg)) lwd <- 7 else lwd <- argg$lwd
  if (!"lend" %in% names(argg)) lend <- 1 else lend <- argg$lend
  if (!"xlim" %in% names(argg)) xlim <- c(0.5,length(x)+0.5) else xlim <- argg$xlim
  if (!"ylim" %in% names(argg)) ylim <- c(0,1) else ylim <- argg$ylim
  if (!"las" %in% names(argg)) las <- 2 else las <- argg$las
  default_cols <- c(grDevices::grey(0.3), grDevices::grey(0.6), 1:7, rep(grDevices::grey(0.9),7))
  names(default_cols) <- c("M-2", "M-1", paste("M", 0:12, sep="+"))
  col <- default_cols[names(x)]
  if ("col" %in% names(argg)) { col <- argg$col }
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par("mar"=opar$mar))
  par(mar=c(ifelse(xlab=="", 3, 5),4,0,0)+0.3)
  plot(as.numeric(x), type="h", axes=FALSE, lwd=lwd, lend=lend, col=col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  graphics::axis(1, at=1:length(x), labels=names(x), las=las)
  graphics::axis(2, las=las)
}
