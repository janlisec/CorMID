#'@title CorMID.
#'@description \code{CorMID} will compute a MID (Mass Isotopomer Distribution)
#'  based on measured ion intensities in GC-APCI-MS.
#'@details Let's assume we measured the ion intensities of all 3 isotopes of an
#'  individual compound containing 2 carbons and observe a vector of {978,22,0}.
#'  We may calculate the enrichment \strong{E} out of this data, i.e. the relative proportion
#'  of 13C vs total carbon which will amount to about 1.1\% (the natural 13C abundance)
#'  under standard conditions.
#'  The equvivalent corMID vector would be {1,0,0}, indicating that the non-labeled
#'  isotopologue (where non-labeled means non-labled above the natural 1.1\%)
#'  is the only component observed.
#'  During a labelling experiment we may change the measurement values in different
#'  ways (either labelling only one carbon or both), which potentially can translate
#'  into similar values for \strong{E} being larger 1.1\%.
#'  The MIDs will provide additional information about the isotopologue fraction
#'  which gave rise to the observed \strong{E}'s (cf. examples). The \emph{r} parameter
#'  indicates an overlay of chemical rearrangements which may occur.
#'@param int Named numeric vector of measured ion intensities of a fragment. Names will give position of values relative to M+H (see details).
#'@param fml Chemical formula of the fragment as string.
#'@param r Either a character vector giving fragments to be considered
#'  OR a named numeric giving relative amounts of fragments
#'  OR NULL (all known fragments will  be estimated)
#'  OR a 2-row matrix giving the lower and upper allowed ratio (see examples).
#'@param penalize Numeric exponent penalizing solutions with low M+H occurence. Formula is 1+3*(1-x)^penalty. NA to omit penalizing.
#'@param mid_fix May provide a numeric vector used as a given MID. Allows to estimate \emph{r} individually.
#'@param trace_steps For testing purposes. Print the results of intermediate steps to console.
#'@param prec Precision of the estimation of MID, set to 1\% as default.
#'@return Estimated percent representation of each isotopologue measured (corMID).
#'@examples
#'# Pyruvic acid 2TMS with 3 biological carbon and 4 ion intensities measured
#'fml <- "C9H20O3Si2"
#'attr(fml,"nmz") <- 7
#'attr(fml,"nbio") <- 3
#'
#'# make up some fake measurement data for this based on assumed 10% labelling at M3
#'mid <- c(0.9,0,0,0.1)
#'r <- unlist(list("M+H"=0.8, "M+"=0.1, "M+H2O-CH4"=0.1))
#'int <- CorMID:::ReconstructMID(mid=mid, r=r, fml=fml)
#'CorMID:::plot_rawMID(int)
#'
#'# full estimation of M and r
#'CorMID(int=int, fml=fml)
#'
#'# get an improved result setting r to the correct values
#'CorMID(int=int, fml=fml, r=r, prec=0.0001)
#'
#'# provoke a wrong estimation using a fixed r
#'CorMID(int=int, fml=fml, r=unlist(list("M+H"=1)))
#'
#'# calculate r if you know the true corMID for a compound
#'r <- attr(CorMID(int=int, fml=fml, mid_fix=c(0.9,0,0,0.1)), "ratio")
#'round(CorMID(int=int, fml=fml, r=r, prec=0.0001),3)
#'
#'# deal with missing intensity values
#'CorMID(int=int[-3], fml=fml)
#'
#'\dontrun{
#'# perform estimation with banded r and observation of optimization steps
#'r <- matrix(c(0.5,1,0,0.5,0,0.5), nrow=2, dimnames=list(NULL,c("M+H","M+","M+H2O-CH4")))
#'CorMID(int=int, fml=fml, r=r, trace=TRUE)
#'
#'#process Gln data from publication
#'utils::data("prep", package = "CorMID")
#'int <- prep[[24]][["int"]][,6]
#'fml <- prep[[24]]$fml
#'CorMID(int=int, fml=fml, trace=TRUE)
#'
#'# check the effect of the penalize paramter on selection of adducts
#'int <- c(1560, 119203, 41927, 16932, 4438)
#'names(int) <- c(-2, 0, 1, 2, 3)
#'fml <- "C19H37NO4Si3"
#'CorMID(int=int, fml=fml, r=NULL, trace=TRUE)
#'CorMID(int=int, fml=fml, r=NULL, trace=TRUE, penalize=7)
#'}
#'@export
CorMID <- function(int=NULL, fml="", r=NULL, penalize=NA, mid_fix=NULL, trace_steps=FALSE, prec = 0.01) {

  # potential parameters
  # specify known fragments; I use this named vector to check in the fitting function where a certain fragment starts,
  # to make the function really flexible this would need to be user definable potentially in a data.frame together with parameter 'r'
  known_frags <- unlist(list("M+H"=0,"M+"=-1,"M-H"=-2,"M+H2O-CH4"=+2))

  # if r is unspecified
  if (is.null(r)) r <- matrix(rep(c(0,1), length(known_frags)), nrow=2, dimnames=list(c("low","high"), names(known_frags)))

  # if r specifies only fragments
  if (is.character(r) && all(r %in% names(known_frags))) r <- matrix(rep(c(0,1),length(r)), nrow=2, dimnames=list(c("low","high"), r))

  # convert r into matrix if numeric vector is provided
  if (is.null(nrow(r)) || nrow(r)==1) {
    r <- r/sum(r) # ensure that sum(r)==1
    r <- matrix(c(r,r), byrow=T, nrow=2, dimnames=list(c("low","high"),names(r)))
  }

  # QC for int and frag
  if (is.null(names(int))) {
    # set frag vector to standard values if not provided by user
    frag <- 0:(length(int)-1)
  } else {
    frag <- as.numeric(gsub("M","",names(int)))
  }
  # ensure that int vector is always complete (fill with zero intensity values if needed)
  if(any(diff(frag)>1)) {
    tmp <- rep(0,diff(range(frag))+1)
    names(tmp) <- frag[1]:frag[length(frag)]
    tmp[as.character(frag)] <- int
    int <- tmp
    frag <- frag[1]:frag[length(frag)]
    attr(fml,"nmz") <- length(int)-1
  }

  # QC for r
  if (!(2 %in% frag) & "M+H2O-CH4" %in% colnames(r)) r <- r[,-which(colnames(r)=="M+H2O-CH4"),drop=FALSE]
  if (!(-2 %in% frag) & "M-H" %in% colnames(r)) r <- r[,-which(colnames(r)=="M-H"),drop=FALSE]
  if (!(-1 %in% frag) & "M+" %in% colnames(r)) r <- r[,-which(colnames(r)=="M+"),drop=FALSE]

  stopifnot(is.matrix(r))
  stopifnot(all(r>=0 & r<=1))

  # QC for fml
  if (is.null(attributes(fml)) && any(frag<0)) {
    cce <- CountChemicalElements(x=fml)
    nC <- cce["C"]-ifelse("Si" %in% names(cce), 3*cce["Si"], 0)
    attr(fml, "nmz") <- nC+abs(min(frag))
    attr(fml, "nbio") <- nC
  }

  # set non-finite intensity values to zero
  if (any(!is.finite(int))) int[!is.finite(int)] <- 0

  # compute theoretical distribution matrix from formula (assuming 1% 13C abundance)
  attr(fml, "nmz") <- length(int)-1
  td <- CalcTheoreticalMDV(fml=fml)

  # QC for mid_fix
  if (!is.null(mid_fix)) {
    if (length(mid_fix)>nrow(td)) mid_fix <- mid_fix[1:nrow(td)]
    if (length(mid_fix)<nrow(td)) warning("length(mid_fix)<nrow(td)")
    names(mid_fix) <- rownames(td)
    mid_fix <- mid_fix/sum(mid_fix)
  }

  # make intensity vector finite (substitute 0 against 1) to improve error calculation
  #int[int==0] <- 1
  # ensure that intensity vector is normalized to sum
  int <- int/sum(int)

  # fitted return value
  out <- FitMID(md=int, td=td, r=r, mid_fix=mid_fix, prec=prec, trace_steps=trace_steps, penalize=penalize)

  # return relative representation of each isotopologue measured (= corrected MIDs in %)
  return(out)
}
