// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// clustering_core
arma::Mat<arma::uword> clustering_core(arma::Mat<unsigned short> m, arma::uword nstates, bool wrap, bool use_8_nb);
RcppExport SEXP _spatialwarnings_clustering_core(SEXP mSEXP, SEXP nstatesSEXP, SEXP wrapSEXP, SEXP use_8_nbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned short> >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nstates(nstatesSEXP);
    Rcpp::traits::input_parameter< bool >::type wrap(wrapSEXP);
    Rcpp::traits::input_parameter< bool >::type use_8_nb(use_8_nbSEXP);
    rcpp_result_gen = Rcpp::wrap(clustering_core(m, nstates, wrap, use_8_nb));
    return rcpp_result_gen;
END_RCPP
}
// pair_counts_internal
arma::Mat<arma::uword> pair_counts_internal(arma::Mat<unsigned short> m, arma::uword nstates, bool wrap, bool use_8_nb);
RcppExport SEXP _spatialwarnings_pair_counts_internal(SEXP mSEXP, SEXP nstatesSEXP, SEXP wrapSEXP, SEXP use_8_nbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<unsigned short> >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nstates(nstatesSEXP);
    Rcpp::traits::input_parameter< bool >::type wrap(wrapSEXP);
    Rcpp::traits::input_parameter< bool >::type use_8_nb(use_8_nbSEXP);
    rcpp_result_gen = Rcpp::wrap(pair_counts_internal(m, nstates, wrap, use_8_nb));
    return rcpp_result_gen;
END_RCPP
}
// coarse_grain_cpp
NumericMatrix coarse_grain_cpp(NumericMatrix mat, int subsize);
RcppExport SEXP _spatialwarnings_coarse_grain_cpp(SEXP matSEXP, SEXP subsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type subsize(subsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(coarse_grain_cpp(mat, subsize));
    return rcpp_result_gen;
END_RCPP
}
// fl_internal
double fl_internal(arma::mat m);
RcppExport SEXP _spatialwarnings_fl_internal(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(fl_internal(m));
    return rcpp_result_gen;
END_RCPP
}
// label_cpp
IntegerMatrix label_cpp(IntegerMatrix mat, IntegerMatrix nbmask, bool wrap);
RcppExport SEXP _spatialwarnings_label_cpp(SEXP matSEXP, SEXP nbmaskSEXP, SEXP wrapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type nbmask(nbmaskSEXP);
    Rcpp::traits::input_parameter< bool >::type wrap(wrapSEXP);
    rcpp_result_gen = Rcpp::wrap(label_cpp(mat, nbmask, wrap));
    return rcpp_result_gen;
END_RCPP
}
// raw_moran
double raw_moran(arma::mat& mat);
RcppExport SEXP _spatialwarnings_raw_moran(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(raw_moran(mat));
    return rcpp_result_gen;
END_RCPP
}
// tplsum
arma::vec tplsum(double expo, double rate, arma::ivec xs, int xmin);
RcppExport SEXP _spatialwarnings_tplsum(SEXP expoSEXP, SEXP rateSEXP, SEXP xsSEXP, SEXP xminSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type expo(expoSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< int >::type xmin(xminSEXP);
    rcpp_result_gen = Rcpp::wrap(tplsum(expo, rate, xs, xmin));
    return rcpp_result_gen;
END_RCPP
}
// tplinfsum
double tplinfsum(double expo, double rate, double xmin, arma::uword maxit, double reltol);
RcppExport SEXP _spatialwarnings_tplinfsum(SEXP expoSEXP, SEXP rateSEXP, SEXP xminSEXP, SEXP maxitSEXP, SEXP reltolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type expo(expoSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    rcpp_result_gen = Rcpp::wrap(tplinfsum(expo, rate, xmin, maxit, reltol));
    return rcpp_result_gen;
END_RCPP
}
// shuffle_matrix
arma::mat shuffle_matrix(arma::mat& mat);
RcppExport SEXP _spatialwarnings_shuffle_matrix(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(shuffle_matrix(mat));
    return rcpp_result_gen;
END_RCPP
}
// shuffle_and_compute
List shuffle_and_compute(arma::mat& mat, Function indic, int nrep);
RcppExport SEXP _spatialwarnings_shuffle_and_compute(SEXP matSEXP, SEXP indicSEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< Function >::type indic(indicSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(shuffle_and_compute(mat, indic, nrep));
    return rcpp_result_gen;
END_RCPP
}
// rspectrum
DataFrame rspectrum(arma::mat mat);
RcppExport SEXP _spatialwarnings_rspectrum(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rspectrum(mat));
    return rcpp_result_gen;
END_RCPP
}
// cpp_skewness
double cpp_skewness(arma::vec X);
RcppExport SEXP _spatialwarnings_cpp_skewness(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_skewness(X));
    return rcpp_result_gen;
END_RCPP
}
// variogram_internal_cpp
NumericMatrix variogram_internal_cpp(NumericMatrix mat, int nmax, int bins, double cutoff);
RcppExport SEXP _spatialwarnings_variogram_internal_cpp(SEXP matSEXP, SEXP nmaxSEXP, SEXP binsSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type nmax(nmaxSEXP);
    Rcpp::traits::input_parameter< int >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(variogram_internal_cpp(mat, nmax, bins, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// sum_all_one_over_k
double sum_all_one_over_k(int from, int to, double expo);
RcppExport SEXP _spatialwarnings_sum_all_one_over_k(SEXP fromSEXP, SEXP toSEXP, SEXP expoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    Rcpp::traits::input_parameter< double >::type expo(expoSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_all_one_over_k(from, to, expo));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatialwarnings_clustering_core", (DL_FUNC) &_spatialwarnings_clustering_core, 4},
    {"_spatialwarnings_pair_counts_internal", (DL_FUNC) &_spatialwarnings_pair_counts_internal, 4},
    {"_spatialwarnings_coarse_grain_cpp", (DL_FUNC) &_spatialwarnings_coarse_grain_cpp, 2},
    {"_spatialwarnings_fl_internal", (DL_FUNC) &_spatialwarnings_fl_internal, 1},
    {"_spatialwarnings_label_cpp", (DL_FUNC) &_spatialwarnings_label_cpp, 3},
    {"_spatialwarnings_raw_moran", (DL_FUNC) &_spatialwarnings_raw_moran, 1},
    {"_spatialwarnings_tplsum", (DL_FUNC) &_spatialwarnings_tplsum, 4},
    {"_spatialwarnings_tplinfsum", (DL_FUNC) &_spatialwarnings_tplinfsum, 5},
    {"_spatialwarnings_shuffle_matrix", (DL_FUNC) &_spatialwarnings_shuffle_matrix, 1},
    {"_spatialwarnings_shuffle_and_compute", (DL_FUNC) &_spatialwarnings_shuffle_and_compute, 3},
    {"_spatialwarnings_rspectrum", (DL_FUNC) &_spatialwarnings_rspectrum, 1},
    {"_spatialwarnings_cpp_skewness", (DL_FUNC) &_spatialwarnings_cpp_skewness, 1},
    {"_spatialwarnings_variogram_internal_cpp", (DL_FUNC) &_spatialwarnings_variogram_internal_cpp, 4},
    {"_spatialwarnings_sum_all_one_over_k", (DL_FUNC) &_spatialwarnings_sum_all_one_over_k, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatialwarnings(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
