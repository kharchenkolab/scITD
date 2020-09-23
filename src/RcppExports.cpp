// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// colMeanVars
Rcpp::DataFrame colMeanVars(SEXP sY, SEXP rowSel, int ncores);
RcppExport SEXP _scITD_colMeanVars(SEXP sYSEXP, SEXP rowSelSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sY(sYSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rowSel(rowSelSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(colMeanVars(sY, rowSel, ncores));
    return rcpp_result_gen;
END_RCPP
}
// get_means
NumericMatrix get_means(SEXP sY, IntegerVector rowSel, IntegerVector numCells);
RcppExport SEXP _scITD_get_means(SEXP sYSEXP, SEXP rowSelSEXP, SEXP numCellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sY(sYSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rowSel(rowSelSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type numCells(numCellsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_means(sY, rowSel, numCells));
    return rcpp_result_gen;
END_RCPP
}
// get_sums
NumericMatrix get_sums(SEXP sY, IntegerVector rowSel);
RcppExport SEXP _scITD_get_sums(SEXP sYSEXP, SEXP rowSelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sY(sYSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rowSel(rowSelSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sums(sY, rowSel));
    return rcpp_result_gen;
END_RCPP
}
// get_vars
NumericVector get_vars(SEXP sY, IntegerVector rowSel, IntegerVector numCells);
RcppExport SEXP _scITD_get_vars(SEXP sYSEXP, SEXP rowSelSEXP, SEXP numCellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sY(sYSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rowSel(rowSelSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type numCells(numCellsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_vars(sY, rowSel, numCells));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scITD_colMeanVars", (DL_FUNC) &_scITD_colMeanVars, 3},
    {"_scITD_get_means", (DL_FUNC) &_scITD_get_means, 3},
    {"_scITD_get_sums", (DL_FUNC) &_scITD_get_sums, 2},
    {"_scITD_get_vars", (DL_FUNC) &_scITD_get_vars, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scITD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}