#' @title Example data as used in function \code{CorMID}.
#' @docType data
#' @format A list of 36 example metabolites.
#' \describe{
#' \item{name}{The compound name.}
#' \item{int}{A numeric matrix providing peak intensities of isotopologues containing biological carbon (rows) over ten samples (columns).}
#' \item{fml}{A character vector of the chemical formula of the compound amended by attributes for the number  of biological carbons 'nbio' and the number of measured mass isotopologues 'nmz'.}
#' }
#' @usage data(prep)
#' @source A flux experiment on SW480/SW620 cell lines by Inna Zaimenko (IZ_Exp05).
"prep"
