#' @title get_idx_mat.
#' @description Calculate an index matrix for isotope distributions.
#' @details tbd.
#' @param n Number of isotopes of this element.
#' @param N Number of atoms of this element in the chemical formula.
#' @return A index matrix together with the frequency vector as an attribute.
#' @keywords internal
#' @noRd
# system.time(
#   lapply(chemforms, function(fml) {
#     x <- cce(fml)[[1]]
#     lapply(names(x), function(ele) {
#       idx <- which(isotopes$element == ele & isotopes$abundance > cutoff)
#       #get_idx_mat(length(idx), x[ele])
#       #get_idx_mat2(length(idx), x[ele])
#       #get_idx_mat3(length(idx), x[ele], abundance_n = isotopes$abundance[idx])
#       precalc_idx[[length(idx)]][[x[ele]]]
#     })
#   })
# )
get_idx_mat <- function(n, N) {
  #x <- expand.grid(lapply(1:n, function(x) {0:N}), KEEP.OUT.ATTRS = FALSE)
  x <- as.matrix(expand.grid(lapply(1:n, function(x) {0:N}), KEEP.OUT.ATTRS = FALSE))
  #dimnames(x) <- NULL
  x <- x[rowSums(x)==N,, drop=FALSE]
  if (n==1) {
    attr(x, "freq") <- 1
  } else {
    attr(x, "freq") <- apply(x, 1, function(x) { choose(N,max(x)) * choose(N-max(x), max(x[-which.max(x)])) })
  }
  #sum(attr(x, "freq")) == n^N
  return(x)
}

# get_idx_mat2 <- function(n, N) {
#   # helper fnc
#   get_comb2 <- function(N) { cbind(seq(N,0), seq(0,N)) }
#   if (n==1) {
#     out <- N
#     attr(out, "freq") <- 1
#   } else if (n==2) {
#     out <- get_comb2(N)
#     attr(out, "freq") <- choose(N,0:N)
#   } else if (n==3) {
#     out <- as.matrix(cbind("0" = rep(N:0, times=1:(N+1)), plyr::ldply(0:N, get_comb2)))
#     attr(out, "freq") <- apply(out, 1, function(x) { choose(N,max(x)) * choose(N-max(x), max(x[-which.max(x)])) })
#   } else {
#     out <- get_idx_mat(n, N)
#   }
#   #dimnames(out) <- NULL
#   return(out)
# }

# get_idx_mat3 <- function(n, N, abundance_n = NULL, cutoff = 0.0001) {
#   if (n==1) {
#     out <- N
#     attr(out, "freq") <- 1
#     return(out)
#   } else {
#     if (is.null(abundance_n)) { N_max <- rep(N, n) } else { N_max <- sapply(1:n, function(i) { max(which(abundance_n[i]^c(1:N)>cutoff)) }) }
#     x <- expand.grid(lapply(1:n, function(x) { 0:N_max[x] }), KEEP.OUT.ATTRS = FALSE)
#     x <- x[rowSums(x)==N,, drop=FALSE]
#     attr(x, "freq") <- apply(x, 1, function(x) { choose(N,max(x)) * choose(N-max(x), max(x[-which.max(x)])) })
#     return(x)
#   }
# }

