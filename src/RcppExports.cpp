// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// interpolateSortedVectorRcpp_
NumericVector interpolateSortedVectorRcpp_(NumericVector vector, int newSize);
RcppExport SEXP blima_interpolateSortedVectorRcpp_(SEXP vectorSEXP, SEXP newSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vector(vectorSEXP);
    Rcpp::traits::input_parameter< int >::type newSize(newSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(interpolateSortedVectorRcpp_(vector, newSize));
    return rcpp_result_gen;
END_RCPP
}