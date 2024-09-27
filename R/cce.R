#' @title cce.
#' @description Count Chemical Elements in a character vector of formulas fast.
#' @details No testing for any chemical alphabet is performed. No brackets,
#'    i.e. [13]C will be removed prior to counting, all elements are expected
#'    to have a trailing number (i.e. '1' can not be omitted).
#' @param x Vector of chemical formulas.
#' @return A list of named numeric vectors with counts for all contained
#'    elements.
#' @examples
#' # count every element
#' cce("C3H7Cl1")
#' cce(c("C3H7Cl1", "C5H12O6", "S10", "C1H2N3O4P5S6"))
#'
#' @keywords internal
#' @noRd
cce <- function(x) {
  e_name <- strsplit(x, "[[:digit:]]+")
  e_numb <- strsplit(x, "[[:alpha:]]+")
  # strsplit(x, "[ABCDEFGHIJKLMNOPQRSTUVWXYZ]")
  return(mapply(function(x, y) {
    structure(as.numeric(x[-1]), names = y)
  }, e_numb, e_name, SIMPLIFY = FALSE))
}
