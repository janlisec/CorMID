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

#' @title A data frame containing isotope information.
#' @docType data
#' @format A data frame of 5 columns for 308 chemical isotopes.
#' \describe{
#' \item{element}{The element name.}
#' \item{isotope}{The isotope name.}
#' \item{mass}{The absolute mass of this isotope in Dalton.}
#' \item{abundance}{The absolute abundance of this isotope.}
#' \item{ratioC}{The C ratio of this element.}
#' }
#' @source Imported from the enviPat package.
"isotopes"

#' @title Pre-calculated index matrices to speed up calculations in function \code{getMID}.
#' @docType data
#' @format A list of 36 example metabolites.
#' @source Calculated using the internal CorMID function `get_idx_mat`.
"precalc_idx"

#' @title The known fragments and their relative position to \code{[M+H]}.
#' @docType data
#' @format A named vector of the currently evaluated fragments.
#' @source Calculated using the internal CorMID function `get_idx_mat`.
"known_frags"
