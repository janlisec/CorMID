#' @title poss_local.
#' @description \code{poss_local} will compute a matrix of possibilities.
#' @details Within the approximation process we need to check various hypotheses
#'  of MID and r combinations. A non-redundant set of posible combinations can
#'  be computed with this \code{poss_local}.
#' @param vec The starting vector sum(vec) should be 1.
#' @param d The maximum allowed deviation for each element in vec.
#' @param prec recision of allowed errors.
#' @param limits A 2-row matrix with lower and upper boundaries for the result vectors.
#' @param ... Passed to function \code{seq}. Either by or length.out (see examples in test-poss_local.R).
#' @return A matrix with rowSums\~1 and within the limits defined by vec and d.
#' @keywords internal
#' @noRd
poss_local <- function(vec=NULL, d=NULL, prec=0.001, limits=NULL, ...) {
  # helper
  fnc_lim <- function(lp=NULL, r=r) {
    flt <- sapply(1:ncol(lp), function(i) {
      lp[,i]>=(r[1,i]-prec) & lp[,i]<=(r[2,i]+prec)
    })
    if (is.null(nrow(flt))) {
      flt <- all(flt)
    } else {
      flt <- apply(flt, 1, all)
    }
    out <- lp[flt,,drop=FALSE]
    return(out)
  }

  if (d==0) {
    lst <- as.list(vec)
  } else {
    lst <- lapply(vec, function(x) {seq(from = max(0, x-d), to = min(1, x+d), ...)})
  }
  tmp <- expand.grid(lst)
  # sum is often not ident(1), hence the solution to compare with precision
  flt <- abs(rowSums(tmp)-1)<=prec
  if (any(flt)) {
    lp <- tmp[flt,,drop=F]
    if (!is.null(limits)) {
      lp <- fnc_lim(lp = lp, r = limits)
    }
  } else {
    vec <- vec/sum(vec)
    lp <- t(unlist(vec))
  }
  rownames(lp) <- 1:nrow(lp)
  # if (nrow(lp)==0) {
  #   lp <- t(as.matrix(vec/sum(vec)))
  # }
  return(as.matrix(lp))
}

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

#' @title verify_suggested.
#' @description Check if packages are available and stop function otherwise.
#' @param pkg Package names to be checked.
#' @return NULL.
#' @keywords internal
#' @noRd
verify_suggested <- function(pkg) {
  # verify that suggested packages are available
  check_pkg <- sapply(pkg, requireNamespace, quietly = TRUE)
  if (!all(check_pkg)) {
    msg <- paste0(
      "The use of this function requires package", ifelse(sum(!check_pkg)>1, "s", ""),
      paste(names(check_pkg)[!check_pkg], collapse=", "),
      ". Please install."
    )
    stop(msg)
  }
  invisible(NULL)
}

#' @title calc_mid_error.
#' @description \code{calc_mid_error} will compute the error of a theoretical and a estimated mid.
#' @param md Normalized measured intensities
#' @param reconstructed_mid A reconstructed MID based on a true MID and a theoretical distribution.
#' @param best_r A named numeric vector of fragment ratios.
#' @return A numeric vector of length(x).
#' @keywords internal
#' @noRd
calc_mid_error <- function(md=NULL, reconstructed_mid=NULL, best_r=NULL) {
  known_frags <- unlist(list("M+H"=0,"M+"=-1,"M-H"=-2,"M+H2O-CH4"=+2))
  length_md <- length(md)
  T0 <- rep(0, max(known_frags))
  frag <- as.numeric(gsub("M","",names(md)))
  out <- rep(0, length_md)
  names(out) <- names(md)
  for (i in 1:length(best_r)) {
    # compute necessary leading and trailing zeros
    L0 <- rep(0, which(frag==known_frags[names(known_frags) %in% names(best_r)[i]])-1)
    # reconstruct the MID for this adduct
    rMID <- reconstructed_mid*unlist(best_r[i])
    # add the MID for this adduct
    out <- out + c(L0, rMID, T0)[1:length_md]
  }
  # normalized and calculate error for this mid
  # $$ the following line led to serious problems and was outcommented on 2022-03-16 by JL
  #if (sum(out)!=1 & sum(out)!=0) out <- out/sum(out)
  mid_err <- sqrt(sum((out-md)^2))
  return(mid_err)
}

#' @title weight_errors.
#' @description \code{weight_errors} will apply a weighting function to a numeric vector of errors.
#' @param rMpH A numeric vector of fragment ratios for the M+H fragment.
#' @param errs A numeric vector of errors.
#' @param penalize A single numeric or NA to omit weighting.
#' @param max_panel Constant to modify weighting function
#' @return A numeric vector of length(x).
#' @keywords internal
#' @noRd
weight_errors <- function(rMpH=NULL, errs=NULL, penalize=NA, max_panel=9) {
  if (is.finite(penalize)) {
    fac <- 1+max_panel*(1-rMpH)^penalize
    return(errs*fac)
  } else {
    return(errs)
  }
}

#' @title select_mid_start_manual.
#' @description \code{select_mid_start_manual} will read in user input interactively
#' @param test_mid Input from `FitMID` to generate messages and perform calculations.
#' @param mid_local Input from `FitMID` to generate messages and perform calculations.
#' @param w_m_errs Input from `FitMID` to generate messages and perform calculations.
#' @return A character vector giving the user selection.
#' @keywords internal
#' @noRd
select_mid_start_manual <- function(test_mid, mid_local, w_m_errs) {
  stopifnot(all(names(test_mid) == names(w_m_errs)))
  stopifnot(all(rownames(mid_local) == names(w_m_errs)))
  cat("\nTop candidates found:\n")
  spacer_col <- data.frame("___"=" | ", check.names = FALSE)
  tmp_print <- cbind(
    round(100*mid_local,2),
    spacer_col,
    t(sapply(test_mid, function(x) { c(formatC(x = x[["r"]], format = "f", digits = 2), spacer_col, "error" = formatC(x[["err"]], format="e", digits=2)) }))
  )
  tmp_print <- cbind(tmp_print, w_m_errs)
  print(utils::head(tmp_print[order(w_m_errs),]))
  # the interactive statement is required to allow a testthat function for this part of FitMID
  if (interactive()) {
    selected_value <- readline(prompt="Type [row_number+enter] to continue stepwise or [enter] without any number to continue to end using best option:")
  } else {
    selected_value <- ""
  }
  if (selected_value=="") {
    # stop calling this function again from 'FitMID' by setting 'trace_steps' to FALSE using an assign statement
    # $$JL$$ assign works as expected but get0() shows that 'trace_steps' is already at FALSE in the outer function (which is not true)
    # assign(x = "trace_steps", value = FALSE, inherits = TRUE, pos = 1)
    # return(mid_local[which.min(w_m_errs),])
    return(NULL)
  } else {
    return(selected_value)
  }
}

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
