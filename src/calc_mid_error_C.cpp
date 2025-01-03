#include <Rcpp.h>
using namespace Rcpp;

//' @description This is an implementation in C to calculate the error between
//' a reconstructed MID and the observed intensity distribution.
//' @keywords internal
//' @noRd

// [[Rcpp::export]]
double calc_mid_error_C(NumericVector md, NumericVector reconstructed_mid, NumericVector best_r, NumericVector r_shift) {
  int n = md.size();
  int m = best_r.size();
  NumericVector out(n);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n-r_shift[j]; ++i) {
      out[i+r_shift[j]] += reconstructed_mid[i] * best_r[j];
    }
  }
  return sqrt(sum(pow(out - md, 2.0)));
}
