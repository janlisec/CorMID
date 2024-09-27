#' @title CorMID.
#' @description \code{CorMID} will compute a MID (Mass Isotopologues
#'     Distribution) based on measured ion intensities in GC-APCI-MS.
#' @details Let's assume we measured the ion intensities of all 3 isotopes of
#'     an individual compound containing 2 carbons and observe a vector
#'     of \{978,22,0\}.
#'     We may calculate the enrichment \strong{E} out of this data, i.e. the
#'     relative proportion of \ifelse{html}{\out{<sup>13</sup>C}}{\eqn{^13^C}}
#'     vs total carbon which will amount to about 1.1\% (the natural
#'     \ifelse{html}{\out{<sup>13</sup>C}}{\eqn{^13^C}} abundance) under
#'     standard conditions.
#'     The equivalent corMID vector would be \{1,0,0\}, indicating that the
#'     non-labeled isotopologue (where non-labeled means non-labeled above the
#'     natural 1.1\%) is the only component observed.
#'     During a labeling experiment we may change the measurement values in
#'     different ways (either labeling only one carbon or both), which
#'     potentially can translate into similar values for \strong{E} being
#'     larger 1.1\%.  The MIDs will provide additional information about the
#'     isotopologue fraction which gave rise to the observed \strong{E}'s (cf.
#'     examples). The \emph{r} parameter  indicates an overlay of chemical
#'     rearrangements which may occur.
#' @param int Named numeric vector of measured ion intensities of a fragment.
#'     Names will give position of values relative to M+H (see details).
#' @param fml Chemical formula of the fragment as string.
#' @param r Either a character vector giving fragments to be considered
#'     OR a named numeric giving relative amounts of fragments
#'     OR NULL (all known fragments will  be estimated)
#'     OR a 2-row matrix giving the lower and upper allowed ratio (see examples).
#' @param penalize Numeric exponent penalizing solutions with low M+H occurrence.
#'     Formula is 1+3*(1-x)^penalty. Set to NA to omit penalizing.
#' @param mid_fix May provide a numeric vector used as a given MID. Allows to
#'     estimate \emph{r} individually.
#' @param trace_steps For testing purposes. Print the results of intermediate
#'     steps to console.
#' @param prec Precision of the estimation of MID, set to 1\% as default.
#' @param algo The algorithm used to estimate the isotopic distribution of a
#'     chemical formula.
#' @return Estimated percent representation of each isotopologue measured (corMID).
#' @references <doi:10.3390/metabo12050408>
#' @examples
#' # make up some fake measurement data for Pyruvic acid 2TMS with 3 biological carbon
#' # assuming 10% labeling at M3 and 2 fragments
#' fml <- "C9H20O3Si2"
#' mid <- c(0.9, 0, 0, 0.1)
#' r <- unlist(list("M+H" = 0.8, "M+H2O-CH4" = 0.2))
#' int <- CorMID::recMID(mid = mid, r = r, fml = fml)
#' plot(int)
#'
#' # full estimation of M and r
#' out <- CorMID::CorMID(int = int, fml = fml)
#' out
#' plot(out)
#'
#' # get an improved result setting r to the correct values
#' CorMID::CorMID(int = int, fml = fml, r = r, prec = 0.0001)
#'
#' # provoke a wrong estimation using a fixed r
#' CorMID::CorMID(int = int, fml = fml, r = unlist(list("M+H" = 1)))
#'
#' # calculate r if you know the true corMID for a compound
#' r <- attr(CorMID::CorMID(int = int, fml = fml, mid_fix = c(0.9, 0, 0, 0.1)), "ratio")
#' round(CorMID::CorMID(int = int, fml = fml, r = r, prec = 0.0001), 3)
#'
#' # deal with missing intensity values
#' CorMID::CorMID(int = int[-3], fml = fml)
#'
#' \donttest{
#' # perform estimation with banded r and observation of optimization steps
#' r <- matrix(c(0.5, 1, 0, 0.5, 0, 0.5), nrow = 2, dimnames = list(NULL, c("M+H", "M+", "M+H2O-CH4")))
#' CorMID::CorMID(int = int, fml = fml, r = r, trace = TRUE)
#'
#' # process Gln data from publication
#' utils::data("prep", package = "CorMID")
#' int <- prep[[24]][["int"]][, 6]
#' fml <- prep[[24]]$fml
#' CorMID::CorMID(int = int, fml = fml, trace = TRUE)
#'
#' # check the effect of the penalize parameter on selection of adducts
#' int <- c(1560, 119203, 41927, 16932, 4438)
#' names(int) <- c(-2, 0, 1, 2, 3)
#' fml <- "C19H37NO4Si3"
#' CorMID::CorMID(int = int, fml = fml, r = NULL, trace = TRUE)
#' CorMID::CorMID(int = int, fml = fml, r = NULL, trace = TRUE, penalize = 7)
#' }
#' @export
CorMID <- function(int = NULL, fml = "", r = NULL, penalize = 7, mid_fix = NULL, trace_steps = FALSE, prec = 0.01, algo = c("CorMID", "Rdisop")) {
  # potential parameters
  # specify known fragments; I use this named vector to check in the fitting function where a certain fragment starts,
  # to make the function really flexible this would need to be user definable potentially in a data.frame together with parameter 'r'
  known_frags <- unlist(list("M+H" = 0, "M+" = -1, "M-H" = -2, "M+H2O-CH4" = +2))
  # this is the maximum number of carbons computable (otherwise matrices get to big)
  lim_nbio <- min(12, attr(fml, "nmz"))

  algo <- match.arg(algo)
  if (algo=="Rdisop") verify_suggested("Rdisop")

  # QC for fml
  # we need to know the biological carbon within fml and specify nmz based on our known frags
  if (is.null(attr(fml, "nbio"))) {
    cce <- CountChemicalElements(x = fml, ele = c("C", "Si"))
    attr(fml, "nbio") <- unname(cce["C"] - 3 * cce["Si"])
  }
  attr(fml, "nbio") <- min(lim_nbio, attr(fml, "nbio"))
  attr(fml, "nmz") <- attr(fml, "nbio") + diff(range(known_frags))
  # compute theoretical distribution matrix from formula (assuming 1% 13C abundance)
  td <- CalcTheoreticalMDV(fml = fml, nbio = attr(fml, "nbio"), nmz = attr(fml, "nmz"), algo = algo)

  # QC for r
  # if r is unspecified
  if (is.null(r)) r <- matrix(rep(c(0, 1), length(known_frags)), nrow = 2, dimnames = list(c("low", "high"), names(known_frags)))
  # if r specifies only fragments
  if (is.character(r) && all(r %in% names(known_frags))) r <- matrix(rep(c(0, 1), length(r)), nrow = 2, dimnames = list(c("low", "high"), r))
  # convert r into matrix if numeric vector is provided
  if (is.null(nrow(r)) || nrow(r) == 1) {
    r <- r / sum(r) # ensure that sum(r)==1
    r <- matrix(c(r, r), byrow = T, nrow = 2, dimnames = list(c("low", "high"), names(r)))
  }

  # QC for int and frag
  frag <- min(known_frags):(max(known_frags) + attr(fml, "nbio"))
  rawMID <- rep(0, length(frag))
  names(rawMID) <- paste0("M", formatC(x = frag, format = "d", flag = "+"))

  # set non-finite intensity values to zero
  if (any(!is.finite(int))) int[!is.finite(int)] <- 0
  # adjust names of int to match rawMID
  if (is.null(names(int))) {
    names(int) <- paste0("M", formatC(x = 0:(length(int) - 1), format = "d", flag = "+"))
  } else {
    tmp <- as.numeric(gsub("M", "", names(int)))
    if (all(is.finite(tmp))) {
      names(int) <- paste0("M", formatC(x = tmp, format = "d", flag = "+"))
    }
  }
  if (!all(names(int) %in% names(rawMID))) stop("rawMID specified without names indicating position relative to [M+H].")
  rawMID[names(int)] <- int

  # QC for mid_fix
  if (!is.null(mid_fix)) {
    if (length(mid_fix) > nrow(td)) mid_fix <- mid_fix[1:nrow(td)]
    if (length(mid_fix) < nrow(td)) warning("length(mid_fix)<nrow(td)")
    names(mid_fix) <- rownames(td)
    mid_fix <- mid_fix / sum(mid_fix)
  }

  # fitted return value
  out <- FitMID(md = rawMID, td = td, r = r, mid_fix = mid_fix, prec = prec, trace_steps = trace_steps, penalize = penalize)

  # set the class to allow dedicated plotting and printing
  class(out) <- "CorMID"

  # return relative representation of each isotopologue measured (= corrected MIDs in %)
  return(out)
}

#' @rdname CorMID
#' @param x Object of class CorMID.
#' @param ... Further plotting parameters.
#' @importFrom graphics axis par
#' @importFrom grDevices grey
#' @export
#' @method plot CorMID
plot.CorMID <- function(x, ...) {
  stopifnot(is.numeric(x))
  stopifnot(all(names(x) %in% paste0("M", 0:12)))
  argg <- c(as.list(environment()), list(...))
  if (!"xlab" %in% names(argg)) xlab <- "" else xlab <- argg$xlab
  if (!"ylab" %in% names(argg)) ylab <- "Relative intensity [%]" else ylab <- argg$ylab
  if (!"lwd" %in% names(argg)) lwd <- 7 else lwd <- argg$lwd
  if (!"lend" %in% names(argg)) lend <- 1 else lend <- argg$lend
  if (!"xlim" %in% names(argg)) xlim <- c(0.5, length(x) + 0.5) else xlim <- argg$xlim
  if (!"ylim" %in% names(argg)) ylim <- c(0, 100) else ylim <- argg$ylim
  if (!"las" %in% names(argg)) las <- 2 else las <- argg$las
  default_cols <- c(1:7, rep(grDevices::grey(0.9), 6))
  names(default_cols) <- paste0("M", 0:12)
  col <- default_cols[names(x)]
  if ("col" %in% names(argg)) {
    col <- argg$col
  }
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par("mar" = opar$mar))
  par(mar = c(ifelse(xlab == "", 3, 5), 4, 0, 0) + 0.3)
  plot(as.numeric(x), type = "h", axes = FALSE, lwd = lwd, lend = lend, col = col, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
  graphics::axis(1, at = 1:length(x), labels = names(x), las = las)
  graphics::axis(2, las = las)
}

#' @rdname CorMID
#' @param x Object of class CorMID.
#' @export
#' @method print CorMID
print.CorMID <- function(x, ...) {
  stopifnot(is.numeric(x))
  # $JL$ I had implemented a colored version of the print to console
  # however, I decided to keep it in standard message style for time being
  cat_or_message <- function(txt) {
    if (interactive() & FALSE) {
      cat(paste0("\033[0;34m", txt, "\033[0m", "\n"))
    } else {
      message(txt)
    }
  }
  cat_or_message("[class] 'CorMID'")
  cat_or_message(paste0("MID [%] (", attr(x, "mid_status"), ")"))
  cat("    ", paste(names(x), collapse = "    "), "\n", sep = "")
  cat(" ", paste(formatC(x, digits = 2, format = "f", width = 5, flag = "0"), collapse = " "), "\n", sep = "")
  cat_or_message(paste0("[attr] 'r' (", attr(x, "ratio_status"), ")"))
  r <- attr(x, "ratio")
  cat("  ", paste(names(r), collapse = "   "), "\n", sep = "")
  cat(sapply(1:length(r), function(i) {
    formatC(r[i], digits = 2, format = "f", width = 2 + nchar(names(r)[i]))
  }), "\n")
  cat_or_message("[attr] 'err'")
  cat(ifelse(is.na(attr(x, "err")), NA, formatC(attr(x, "err"), format = "g")), "\n\n")
}
