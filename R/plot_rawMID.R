#'@title plot_rawMID.
#'@description \code{plot_rawMID} will plot a raw mid with some color standards.
#'@param x A raw mid vector.
#'@param ... Used to pass arguments for plot to overwrite some default settings.
#'@return NULL (plot is generated).
#'@keywords internal
#'@importFrom graphics axis par
#'@importFrom grDevices grey
#'@examples
#'fml <- "C9H20O3Si2"
#'attr(fml,"nmz") <- 5
#'attr(fml,"nbio") <- 3
#'# make up some fake measurement data for this based on assumed 10% labelling at M3
#'mid <- c(0.9,0,0,0.1)
#'r <- unlist(list("M+H"=0.8, "M+"=0.1, "M+H2O-CH4"=0.1))
#'rMID <- CorMID:::ReconstructMID(mid=mid, r=r, fml=fml)
#'CorMID:::plot_rawMID(x = rMID)
#'CorMID:::plot_rawMID(x = rMID, ylim=c(0,max(rMID)))
#'CorMID:::plot_rawMID(x = rMID, xlim=c(-2,12), ylim=NULL, col=2, lwd=12, las=2, xlab="label")
plot_rawMID <- function(x = NULL, ...) {
  stopifnot(is.numeric(x))
  stopifnot(is.vector(x))
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
  if (!"las" %in% names(argg)) las <- 1 else las <- argg$las
  default_cols <- c(grDevices::grey(0.3),grDevices::grey(0.6),1:7,rep(grDevices::grey(0.9),7))
  names(default_cols) <- c("M-2", "M-1", paste("M", 0:12, sep="+"))
  col <- default_cols[names(x)]
  if ("col" %in% names(argg)) { col <- argg$col }
  opar <- graphics::par(no.readonly = TRUE)
  par(mar=c(ifelse(xlab=="", 3, 5),4,0,0)+0.3)
  plot(x, type="h", axes=FALSE, lwd=lwd, lend=lend, col=col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  graphics::axis(1, at=1:length(x), labels=names(x), las=las)
  graphics::axis(2, las=las)
  graphics::par("mar"=opar$mar)
  invisible(NULL)
}
