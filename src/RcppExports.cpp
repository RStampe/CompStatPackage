// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// perform_em_cpp
List perform_em_cpp(NumericVector x, List listOfInitialPar, int maxIteration, double epsilon);
RcppExport SEXP _CompStatPackage_perform_em_cpp(SEXP xSEXP, SEXP listOfInitialParSEXP, SEXP maxIterationSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type listOfInitialPar(listOfInitialParSEXP);
    Rcpp::traits::input_parameter< int >::type maxIteration(maxIterationSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(perform_em_cpp(x, listOfInitialPar, maxIteration, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// kernel_cpp_epanechnikov
List kernel_cpp_epanechnikov(const NumericVector& grid_points, const NumericVector& vec_obs, double bandwidth);
RcppExport SEXP _CompStatPackage_kernel_cpp_epanechnikov(SEXP grid_pointsSEXP, SEXP vec_obsSEXP, SEXP bandwidthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type grid_points(grid_pointsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type vec_obs(vec_obsSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_cpp_epanechnikov(grid_points, vec_obs, bandwidth));
    return rcpp_result_gen;
END_RCPP
}
// kernel_cpp_gaussian
List kernel_cpp_gaussian(const NumericVector& grid_points, const NumericVector& vec_obs, double bandwidth);
RcppExport SEXP _CompStatPackage_kernel_cpp_gaussian(SEXP grid_pointsSEXP, SEXP vec_obsSEXP, SEXP bandwidthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type grid_points(grid_pointsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type vec_obs(vec_obsSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_cpp_gaussian(grid_points, vec_obs, bandwidth));
    return rcpp_result_gen;
END_RCPP
}
// kernel_cpp_epanechnikov_binning
List kernel_cpp_epanechnikov_binning(const NumericVector& grid_points, const NumericVector& vec_obs, double bandwidth);
RcppExport SEXP _CompStatPackage_kernel_cpp_epanechnikov_binning(SEXP grid_pointsSEXP, SEXP vec_obsSEXP, SEXP bandwidthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type grid_points(grid_pointsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type vec_obs(vec_obsSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_cpp_epanechnikov_binning(grid_points, vec_obs, bandwidth));
    return rcpp_result_gen;
END_RCPP
}
// kernel_cpp_gaussian_binning
List kernel_cpp_gaussian_binning(const NumericVector& grid_points, const NumericVector& vec_obs, double bandwidth);
RcppExport SEXP _CompStatPackage_kernel_cpp_gaussian_binning(SEXP grid_pointsSEXP, SEXP vec_obsSEXP, SEXP bandwidthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type grid_points(grid_pointsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type vec_obs(vec_obsSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    rcpp_result_gen = Rcpp::wrap(kernel_cpp_gaussian_binning(grid_points, vec_obs, bandwidth));
    return rcpp_result_gen;
END_RCPP
}
// set_data_vectors
void set_data_vectors(NumericVector x_in, NumericVector z_in);
RcppExport SEXP _CompStatPackage_set_data_vectors(SEXP x_inSEXP, SEXP z_inSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_in(x_inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z_in(z_inSEXP);
    set_data_vectors(x_in, z_in);
    return R_NilValue;
END_RCPP
}
// log_target_density_cpp
NumericVector log_target_density_cpp(NumericVector y);
RcppExport SEXP _CompStatPackage_log_target_density_cpp(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(log_target_density_cpp(y));
    return rcpp_result_gen;
END_RCPP
}
// d_log_target_density_cpp
NumericVector d_log_target_density_cpp(NumericVector y);
RcppExport SEXP _CompStatPackage_d_log_target_density_cpp(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(d_log_target_density_cpp(y));
    return rcpp_result_gen;
END_RCPP
}
// get_accepts_gaussian_cpp
LogicalVector get_accepts_gaussian_cpp(double mu, double sd, double log_alpha_prime, NumericVector proposals, NumericVector u_samples);
RcppExport SEXP _CompStatPackage_get_accepts_gaussian_cpp(SEXP muSEXP, SEXP sdSEXP, SEXP log_alpha_primeSEXP, SEXP proposalsSEXP, SEXP u_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_alpha_prime(log_alpha_primeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proposals(proposalsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u_samples(u_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_accepts_gaussian_cpp(mu, sd, log_alpha_prime, proposals, u_samples));
    return rcpp_result_gen;
END_RCPP
}
// create_sampler_function_cpp
NumericVector create_sampler_function_cpp(int n, double c, const NumericVector& Q, const NumericVector& a, const NumericVector& b, const NumericVector& z);
RcppExport SEXP _CompStatPackage_create_sampler_function_cpp(SEXP nSEXP, SEXP cSEXP, SEXP QSEXP, SEXP aSEXP, SEXP bSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(create_sampler_function_cpp(n, c, Q, a, b, z));
    return rcpp_result_gen;
END_RCPP
}
// create_log_envelope_density_cpp
NumericVector create_log_envelope_density_cpp(const NumericVector& x, const NumericVector& z, const NumericVector& a, const NumericVector& b);
RcppExport SEXP _CompStatPackage_create_log_envelope_density_cpp(SEXP xSEXP, SEXP zSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(create_log_envelope_density_cpp(x, z, a, b));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CompStatPackage_perform_em_cpp", (DL_FUNC) &_CompStatPackage_perform_em_cpp, 4},
    {"_CompStatPackage_kernel_cpp_epanechnikov", (DL_FUNC) &_CompStatPackage_kernel_cpp_epanechnikov, 3},
    {"_CompStatPackage_kernel_cpp_gaussian", (DL_FUNC) &_CompStatPackage_kernel_cpp_gaussian, 3},
    {"_CompStatPackage_kernel_cpp_epanechnikov_binning", (DL_FUNC) &_CompStatPackage_kernel_cpp_epanechnikov_binning, 3},
    {"_CompStatPackage_kernel_cpp_gaussian_binning", (DL_FUNC) &_CompStatPackage_kernel_cpp_gaussian_binning, 3},
    {"_CompStatPackage_set_data_vectors", (DL_FUNC) &_CompStatPackage_set_data_vectors, 2},
    {"_CompStatPackage_log_target_density_cpp", (DL_FUNC) &_CompStatPackage_log_target_density_cpp, 1},
    {"_CompStatPackage_d_log_target_density_cpp", (DL_FUNC) &_CompStatPackage_d_log_target_density_cpp, 1},
    {"_CompStatPackage_get_accepts_gaussian_cpp", (DL_FUNC) &_CompStatPackage_get_accepts_gaussian_cpp, 5},
    {"_CompStatPackage_create_sampler_function_cpp", (DL_FUNC) &_CompStatPackage_create_sampler_function_cpp, 6},
    {"_CompStatPackage_create_log_envelope_density_cpp", (DL_FUNC) &_CompStatPackage_create_log_envelope_density_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_CompStatPackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
