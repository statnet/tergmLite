#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector shiftVec(NumericVector elv, NumericVector vid) {

  int n = elv.size();
  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    out[i] = elv[i] - sum(elv[i] > vid);
  }

  return(out);

}
