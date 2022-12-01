#'@title CalcTheoreticalMDV.
#'@description \code{CalcTheoreticalMDV} will compute the Mass Distribution Vectors
#'  of isotopologues as it is used for correction matrix in \link{CorMID} computations.
#'@details \code{CalcTheoreticalMDV} basically is a convenience function using Rdisop
#'  to generate the isotopologue distribution at natural abundance of 13C for a given
#'  formula.
#'  It will break this down into a matrix where the components of the MID constitute
#'  the rows and the expected relative ion intensities are within the columns.
#'  The number of exported ion intensities and MID components can be limited
#'  if numeric values for "nmz" and/or "nbio" are provided as parameters.
#'@param fml The chemical formula of the compound.
#'@param nbio Provide the number of biological carbon within fml explicitly.
#'@param nmz Provide the number of measured isotopes of fml explicitly.
#'@return A matrix of theoretical mass distribution vectors.
#'@importFrom Rdisop getMolecule
#'@examples
#'# standard distribution matrix
#'fml <- "C5H6Si1"
#'CalcTheoreticalMDV(fml=fml)
#'
#'# extend to more columns (number of measured ions) if required
#'CalcTheoreticalMDV(fml=fml, nmz=4)
#'
#'# limit to a smaller number of biological carbon (i.e. if compounds are silylated)
#'CalcTheoreticalMDV(fml=fml, nmz=4, nbio=2)
#'@export
#CalcTheoreticalMDV <- function(fml=NULL, nbio=NULL, nmz=NULL, algo=c("Rdisop", "IsoSpec")[2]) {
CalcTheoreticalMDV <- function(fml=NULL, nbio=NULL, nmz=NULL) {
  # establish number of ions measured (estimate from formula if not provided)
  cce <- CountChemicalElements(x=fml)
  if (is.null(nmz)) {
    maxisotopes <- cce["C"]-ifelse("Si" %in% names(cce), 3*cce["Si"], 0)
  } else {
    maxisotopes <- nmz
  }
  # establish number of biological carbons (estimate from formula if not provided)
  if (is.null(nbio) || nbio>maxisotopes) {
    maxmid <- maxisotopes
  } else {
    maxmid <- nbio
  }
  # get base distribution and convert into matrix including normalization
  td <- matrix(NA, ncol=maxisotopes+1, nrow=maxmid+1, dimnames=list(paste0("M",0:maxmid),paste0("M+",0:maxisotopes)))
  # if (algo == "Rdisop") {
    for (i in 1:nrow(td)) {
      fml2 <- paste0(names(cce), cce, collapse="")
      bd <- Rdisop::getMolecule(fml2, maxisotopes = maxisotopes+1-(i-1))$isotopes[[1]][2,]
      td[i,] <- c(rep(0, i-1), bd)
      cce["C"] <- cce["C"]-1
    }
  # } else {
  #   # use IsoSpec algo
  #   for (i in 1:nrow(td)) {
  #     # keep precision at 10^-9 to account for mass spec intensity range (probably 10^-6 would be sufficient)
  #     tmp <- IsoSpecR::IsoSpecify(molecule = cce, stopCondition = .99999999)
  #     # combine isotopic fine structure into a isotope vector
  #     split_var <- round(tmp[,1]-min(tmp[,1]))
  #     n <- maxisotopes+1-(i-1)
  #     bd <- tapply(tmp[,2], split_var, sum)
  #     #bd <- tapply(tmp[,2], split_var, max)
  #     if (length(bd)<n) bd <- c(bd, rep(0, n-length(bd)))
  #     if (length(bd)>n) bd <- bd[1:n]
  #     td[i,] <- c(rep(0, i-1), bd)
  #     cce["C"] <- cce["C"]-1
  #   }
  # }
  td <- t(apply(td, 1, function(x) {x/sum(x)}))
  # return matrix
  return(td)
}
