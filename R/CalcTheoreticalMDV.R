#' @title CalcTheoreticalMDV.
#' @description \code{CalcTheoreticalMDV} will compute the Mass Distribution Vectors
#'  of isotopologues as it is used for correction matrix in \link{CorMID} computations.
#' @details \code{CalcTheoreticalMDV} basically is a convenience function using Rdisop
#'  to generate the isotopologue distribution at natural abundance of \eqn{^{13}C}{13C}
#'  for a given formula.
#'  It will break this down into a matrix where the components of the MID constitute
#'  the rows and the expected relative ion intensities are within the columns.
#'  The number of exported ion intensities and MID components can be limited
#'  if numeric values for \code{nmz} and/or \code{nbio} are provided as parameters.
#' @param fml The chemical formula of the compound.
#' @param nbio Provide the number of biological carbon within \code{fml} explicitly.
#' @param nmz Provide the number of measured isotopes of \code{fml} explicitly.
#' @param algo algo.
#' @return A matrix of theoretical mass distribution vectors.
#' @examples
#' # standard distribution matrix
#' fml <- "C5H6Si1"
#' CalcTheoreticalMDV(fml = fml)
#'
#' # extend to more columns (number of measured ions) if required
#' CalcTheoreticalMDV(fml = fml, nmz = 4)
#'
#' # limit to a smaller number of biological carbon (i.e. if compounds are silylated)
#' CalcTheoreticalMDV(fml = fml, nmz = 4, nbio = 2)
#' CalcTheoreticalMDV(fml = fml, nmz = 4, nbio = 2, algo="Rdisop")
#'
#' # from Vignette
#' fml <- "C21Si5"
#' round(CalcTheoreticalMDV(fml = fml, nbio = 21, nmz = 21)[-(5:19), -(5:19)], 4)
#' @export
#CalcTheoreticalMDV <- function(fml=NULL, nbio=NULL, nmz=NULL, algo=c("Rdisop", "IsoSpec")[2]) {
CalcTheoreticalMDV <- function(fml = NULL, nbio = NULL, nmz = NULL, algo = c("CorMID", "Rdisop")) {
  algo <- match.arg(algo)
  if (algo=="Rdisop") verify_suggested("Rdisop")
  # establish number of ions measured (estimate from formula if not provided)
  nce <- CountChemicalElements(x=fml)
  if (is.null(nmz)) {
    maxisotopes <- nce["C"]-ifelse("Si" %in% names(nce), 3*nce["Si"], 0)
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
  #browser()
  if (algo == "Rdisop") {
    for (i in 1:nrow(td)) {
      bd <- Rdisop::getMolecule(paste0(names(nce), nce, collapse=""), maxisotopes = maxisotopes+1-(i-1))$isotopes[[1]][2,]
      td[i,] <- c(rep(0, i-1), bd)
      if ("C" %in% names(nce) && nce["C"]>=2) nce["C"] <- nce["C"]-1 else nce <- nce[names(nce)!="C"]
    }
  # } else if (algo == "IsoSpec") {
  #   # use IsoSpec algo
  #   for (i in 1:nrow(td)) {
  #     # keep precision at 10^-9 to account for mass spec intensity range (probably 10^-6 would be sufficient)
  #     tmp <- IsoSpecR::IsoSpecify(molecule = nce, stopCondition = .99999999)
  #     # combine isotopic fine structure into a isotope vector
  #     split_var <- round(tmp[,1]-min(tmp[,1]))
  #     n <- maxisotopes+1-(i-1)
  #     bd <- tapply(tmp[,2], split_var, sum)
  #     if (length(bd)<n) bd <- c(bd, rep(0, n-length(bd)))
  #     if (length(bd)>n) bd <- bd[1:n]
  #     td[i,] <- c(rep(0, i-1), bd)
  #     if ("C" %in% names(nce) && nce["C"]>=2) nce["C"] <- nce["C"]-1 else nce <- nce[names(nce)!="C"]
  #   }
  } else if (algo == "CorMID") {
    for (i in 1:nrow(td)) {
      bd <- getMID(fml = paste0(names(nce), nce, collapse=""), prec = 8)[,2]
      n <- maxisotopes+1-(i-1)
      if (length(bd)<n) bd <- c(bd, rep(0, n-length(bd)))
      if (length(bd)>n) bd <- bd[1:n]
      td[i,] <- c(rep(0, i-1), bd)
      if ("C" %in% names(nce) && nce["C"]>=2) nce["C"] <- nce["C"]-1 else nce <- nce[names(nce)!="C"]
    }
  }
  td <- t(apply(td, 1, function(x) {x/sum(x)}))
  # return matrix
  return(td)
}
