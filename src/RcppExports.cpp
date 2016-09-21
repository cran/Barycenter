// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Subgradient
arma::vec Subgradient(arma::vec a, arma::vec b, arma::mat M, double lambda, int maxIter, double tolerance);
RcppExport SEXP Barycenter_Subgradient(SEXP aSEXP, SEXP bSEXP, SEXP MSEXP, SEXP lambdaSEXP, SEXP maxIterSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(Subgradient(a, b, M, lambda, maxIter, tolerance));
    return rcpp_result_gen;
END_RCPP
}
