#'@title CalcTheoreticalMDV.
#'@description \code{CalcTheoreticalMDV} will compute the Mass Distribution Vectors
#'  of isotopologues as it is used for correction matrix in \link{CorMID} computations.
#'@details \code{CalcTheoreticalMDV} basically is a convenience function using Rdisop
#'  to generate the isotopologue distribution at natural abundance of 13C for a given
#'  formula.
#'  It will break this down into a matrix where the components of the MID constitute
#'  the rows and the expected relative ion intensities are within the columns.
#'  The number of exported ion intensities and MID components can be limited
#'  if numeric values for "nmz" and/or "nbio" provided as attributes with the formula.
#'@param fml The chemical formula of the compound.
#'@return A matrix of theoretical mass distribution vectors.
#'@importFrom Rdisop getMolecule
#'@examples
#'# standard distribution matrix
#'fml <- "C5H6Si1"
#'CalcTheoreticalMDV(fml=fml)
#'
#'# extend to more columns (number of measured ions) if required
#'attr(fml,"nmz") <- 4
#'CalcTheoreticalMDV(fml=fml)
#'
#'# limit to a smaller number of biological carbon (i.e. if compounds are silylated)
#'attr(fml,"nbio") <- 2
#'CalcTheoreticalMDV(fml=fml)
#'@export
CalcTheoreticalMDV <- function(fml=NULL) {
  # establish number of ions measured (estimate from formula if not provided)
  cce <- CountChemicalElements(x=fml)
  if (is.null(attr(fml, "nmz"))) {
    maxisotopes <- cce["C"]-ifelse("Si" %in% names(cce), 3*cce["Si"], 0)
  } else {
    maxisotopes <- attr(fml, "nmz")
  }
  # establish number of biological carbons (estimate from formula if not provided)
  if (is.null(attr(fml, "nbio")) || attr(fml, "nbio")>maxisotopes) {
    maxmid <- maxisotopes
  } else {
    maxmid <- attr(fml, "nbio")
  }
  # get base distribution and convert into matrix including normalization
  td <- matrix(NA, ncol=maxisotopes+1, nrow=maxmid+1, dimnames=list(paste0("M",0:maxmid),paste0("M+",0:maxisotopes)))
  for (i in 1:nrow(td)) {
    tmp <- paste0(names(cce),cce,collapse="")
    bd <- Rdisop::getMolecule(tmp, maxisotopes = maxisotopes+1-(i-1))$isotopes[[1]][2,]
    td[i,] <- c(rep(0,i-1), bd)
    cce["C"] <- cce["C"]-1
  }
  td <- t(apply(td, 1, function(x) {x/sum(x)}))
  # returrn matrix
  return(td)
}
