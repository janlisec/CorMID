#' @title getMID.
#' @description \code{getMID} will determine the measurable isotopic
#'    spectrum for a chemical formula.
#' @details The computation yields similar results that would be obtained by
#'    packages `Rdisop` or `enviPat` but is completely in R (no C++ dependencies).
#'    However, it is approx. 7-fold slower than `Rdisop`. Where processing speed
#'    is of importance, please use the `algo` parameter of the `CorMID` function.
#' @param fml Chemical formula.
#' @param resolution Currently fixed to 20000 (might be made changeable in the future).
#' @param cutoff Discard peaks below this threshold (relative to highest peak).
#' @param isotopes Specify explicitly or keep NULL to use internally provided list.
#' @param prec Rounding precision of returned mz and int values.
#' @param step Can be used to return intermediate results (might be deprecated in the future).
#' @return A two column matrix for mz and int values of the calculated spectrum.
#' @examples
#' fml <- "C3H7Cl1"
#' getMID(fml)
#' \dontrun{
#' bench::mark(
#'  CorMID = dim(getMID(fml, prec=5)),
#'  Rdisop = dim(round(t(Rdisop::getMolecule(fml)$isotopes[[1]])[1:4,],5))
#' )
#' }
#'
#' \dontrun{
#' data(chemforms, package = "enviPat")
#' chemforms <- chemforms[-grep("[[]", chemforms)]
#' bench::mark(
#'  CorMID = length(lapply(chemforms, getMID)),
#'  Rdisop = length(lapply(chemforms, Rdisop::getMolecule))
#' )
#' }
#' @export
getMID <- function(fml, resolution = 20000, cutoff = 0.0001, isotopes = NULL, prec = 4, step = 0) {
  #print(fml)
  if (is.null(isotopes)) isotopes <- CorMID::isotopes
  if (fml=="") return(data.frame("mz"=12, "int"=1))
  #browser()
  precalc_idx <- CorMID::precalc_idx
  # count elements
  fml_cce <- cce(fml)[[1]]

  # get element individual distributions and limit to cutoff
  ele_mzint <- lapply(names(fml_cce), function(ele) {
    idx <- which(isotopes$element == ele & isotopes$abundance > cutoff)
    # $$JL$$ several on-the-fly versions to compute the index matrix exist, however pre-calculation is still the fastest
    #idx_mat <- get_idx_mat(n, N)
    #idx_mat <- get_idx_mat2(n, N)
    #idx_mat <- get_idx_mat3(n, N, abundance_n = isotopes$abundance[idx])
    idx_mat <- precalc_idx[[length(idx)]][[fml_cce[ele]]]
    out <- data.frame(
      "mz" = colSums(t(idx_mat) * isotopes$mass[idx]),
      "int" = apply(isotopes$abundance[idx]^t(idx_mat), 2, prod) * attr(idx_mat, "freq")
    )
    return(out[out[, 2] > cutoff, ])
  })
  if (step==1) return(ele_mzint)

  # combine elements
  if (length(ele_mzint)==1) {
    df_jl <- ele_mzint[[1]]
  } else {
    idx <- expand.grid(lapply(ele_mzint, function(x) { 1:nrow(x) }))
    df_jl <- data.frame(
      "mz" = rowSums(sapply(1:length(ele_mzint), function(i) { ele_mzint[[i]][idx[, i], 1] })),
      "int" = apply(sapply(1:length(ele_mzint), function(i) { ele_mzint[[i]][idx[, i], 2] }), 1, prod)
    )
  }
  if (step==2) return(df_jl)

  # reduce the peaks according to approx. 20k resolution
  df_jl <- t(sapply(split(df_jl, round(df_jl[, 1]-min(df_jl[, 1]))), function(x) {
    c(stats::weighted.mean(x = x[, 1], w = x[, 2]), sum(x[, 2]))
  }, USE.NAMES = FALSE))
  dimnames(df_jl) <- list(NULL, c("mz","int"))
  if (step==3) return(df_jl)

  #df_jl[, 2] <- df_jl[, 2] / max(df_jl[, 2])
  # round mz and int as specified by prec'
  if (is.finite(prec[1])) df_jl <- round(df_jl, prec[1])
  return(df_jl[df_jl[, 2] > cutoff, , drop=FALSE])
}

# getMID2 <- function(fml, resolution = 20000, cutoff = 0.0001, isotopes = NULL, prec = 4, step = 0) {
#   if (is.null(isotopes)) isotopes <- CorMID::isotopes
#   isotopes[isotopes$element=="D","isotope"] <- "D1"
#   fml_cce <- cce(fml)[[1]]
#   isos <- isotopes[isotopes$element %in% names(fml_cce) & isotopes$abundance > cutoff,]
#   iso_freq <- table(isos$element)
#   isos$n <- fml_cce[isos$element]
#   isos$max <- floor(log10(cutoff)/log10(isos$abundance))
#   isos$max[!is.finite(isos$max)] <- 1
#   isos$min <- 0
#   isos$min[isos$element %in% names(which(iso_freq==1))] <- 1
#
#   if (step==1) return(isos)
#
#   rownames(isos) <- isos$isotope
#
#   # $$JL$$ the expand.grid step is really slooooow for some formulas
#   tmp <- expand.grid(sapply(isos$isotope, function(x) { isos[x,"min"]:min(isos[x,c("n","max")]) }))
#   if (step==2) return(tmp)
#   for (e in names(iso_freq[iso_freq>=2])[order(fml_cce[names(iso_freq[iso_freq>=2])], decreasing = T)]) {
#     ks <- which(e == isos$element)
#     tmp <- tmp[rowSums(tmp[,ks])==fml_cce[e],]
#   }
#
#
#   s <- data.frame(
#     colSums(isos$mass * t(tmp)),
#     apply(isos$abundance^t(tmp), 2, prod)
#   )
#
#   if (step==3) return(s)
#
#   df_jl <- plyr::ldply(split(s, round(s[, 1]-min(s[, 1]))), function(x) {
#     data.frame(
#       "mz" = stats::weighted.mean(x = x[, 1], w = x[, 2]),
#       "int" = sum(x[, 2])
#     )
#   }, .id = NULL)
#
#   if (step==4) return(df_jl)
#
#   if (is.finite(prec[1])) df_jl <- round(df_jl, prec[1])
#   return(df_jl[df_jl[, 2] > cutoff, ])
# }
#
# # bench::mark(
# #   "a" = rowSums(tmp[,1:2])==7,
# #   "b" = (tmp[,1]+tmp[,2])==7
# # )
