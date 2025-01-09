#include <Rcpp.h>
using namespace Rcpp;

//' @title poss_local_C.
//' @description \code{poss_local_C} will compute a matrix of possibilities.
//' @details Within the approximation process we need to check various hypotheses
//'    of MID and r combinations. A non-redundant set of possible combinations can
//'    be computed with this \code{poss_local_C}.
//' @param vec The starting vector sum(vec) should be 1.
//' @param d The maximum allowed deviation for each element in vec.
//' @param prec Precision of allowed errors.
//' @param limits A 2-row matrix with lower and upper boundaries for the result vectors.
//' @param by Passed as parameter to function \code{seq}.
//' @param length Passed as parameter to function \code{seq}. Not evaluated if by is specified.
//' @return A matrix with rowSums\~1 and within the limits defined by vec, d and limits.
//' @keywords internal
//' @noRd

// helper function allowing to efficiently apply limits to a poss_local matrix
NumericMatrix fnc_lim(NumericMatrix lp, NumericMatrix limits)
{
  int nrow = lp.nrow();
  int ncol = lp.ncol();
  NumericVector out;
  bool check_row;
  for (int i = 0; i < nrow; i++) {
    check_row = true;
    for (int j = 0; j < ncol; j++) {
      check_row = check_row && lp(i,j)>=limits(0,j) && lp(i,j)<=limits(1,j);
    }
    if (check_row) {
      for (int j = 0; j < ncol; j++) {
        out.push_back(lp(i,j));
      }
    }
  }
  out.attr("dim") = Dimension(ncol, out.length()/ncol);
  NumericMatrix m = as<NumericMatrix>(out);
  return(transpose(m));
}

NumericVector seqC(double from, double to, double by = 1.0, int length = 0, double min_val = 0, double max_val = 1.0) {
  if (length>0) {
    by = std::abs(to - from)/(length-1);
  } else {
    length = 1+round(std::abs(to - from)/by);
  }
  std::vector<double> out;
  out.reserve(length);
  if (to < from) {
    by = -1 * by;
  }
  for (int i = 0; i < length; i++) {
    if ((from >= min_val) & (from <= max_val)) {
      out.push_back(from);
    }
    from += by;
  }
  return(wrap(out));
}

NumericMatrix rowSums_flt(NumericMatrix x, double prec = 0.001, double exp_sum = 1) {
  int nrow = x.nrow(), ncol = x.ncol();
  std::vector<double> out;
  out.reserve(nrow*ncol);
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    if (std::abs(total - exp_sum)<prec) {
      for (int j = 0; j < ncol; j++) {
        out.push_back(x(i, j));
      }
    }
  }
  NumericVector tmp = wrap(out);
  tmp.attr("dim") = Dimension(ncol, tmp.length()/ncol);
  NumericMatrix out_mat = as<NumericMatrix>(tmp);
  return(transpose(out_mat));
}

NumericMatrix expGrid(List lst) {
  int i; int j; int k; int l;
  int n = lst.length();
  std::vector<int>  m(n);
  std::vector<int> inner_rep(n, 1);
  std::vector<int> outer_rep(n, 1);
  NumericVector tmp;
  std::vector<double> out;
  for (i = 0; i < n; i++) {
    tmp = lst[i];
    m[i] = tmp.length();
  }
  for (i = 1; i < n; i++) {
    inner_rep[i] *= inner_rep[i-1] * m[i-1];
  }
  out.reserve(m[n-1]*inner_rep[n-1]);
  for (i = n-1; i-- > 0;) {
    outer_rep[i] *= outer_rep[i+1] * m[i+1];
  }
  for (i = 0; i < n; i++) {
    tmp = lst[i];
    for (k = 0; k < outer_rep[i]; k++) {
      for (j = 0; j < m[i]; j++) {
        for (l = 0; l < inner_rep[i]; l++) {
          out.push_back(tmp[j]);
        }
      }
    }
  }
  tmp = wrap(out);
  tmp.attr("dim") = Dimension(tmp.length()/n, n);
  NumericMatrix out_mat = as<NumericMatrix>(tmp);
  return(out_mat);
}


// [[Rcpp::export]]
NumericMatrix poss_local_C(NumericVector vec, double d, double prec = 0.001, Nullable<NumericMatrix> limits = R_NilValue, double by = 1.0, int length = 0) {
  int i;
  int n = vec.size();
  CharacterVector cn(n, "Var");
  NumericVector rn;
  if (vec.hasAttribute("names")) {
    cn = vec.names();
  } else {
    for (i = 0; i < n; i++) {
      cn[i] += i;
    }
  }
  NumericMatrix out;
  // establish matrix of possible combinations of r
  List lst(n);
  for (i = 0; i < n; i++) {
    lst[i] = seqC(vec[i]-d, vec[i]+d, by, length);
  }
  out = expGrid(lst);
  // filter for combinations using condition RowSums == 1
  out = rowSums_flt(out);
  // filter using an additional parameter specifying ranges (limits) for each column
  if (limits.isNotNull()) {
    NumericMatrix limits_(limits);
    if (limits_.ncol() == out.ncol()) {
      out = fnc_lim(out, limits_);
    }
  }
  // set matrix rownames
  n = out.nrow();
  rn = seq_len(n);
  rownames(out) = rn;
  colnames(out) = cn;
  return(out);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
poss_local_C(vec=c(0.5,0.25,0.25), d=0.25)
*/
