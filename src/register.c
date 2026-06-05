#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Forward declarations of our Zig functions
extern SEXP zd_euclidean_dist_self_R(SEXP x_sexp, SEXP num_threads_sexp, SEXP has_na_sexp);
extern SEXP zd_euclidean_dist_cross_R(SEXP x_sexp, SEXP y_sexp, SEXP num_threads_sexp, SEXP has_na_sexp);
extern SEXP zd_gower_dist_self_R(SEXP x_num, SEXP x_cat, SEXP ranges, SEXP w_num, SEXP w_cat, SEXP num_threads_sexp, SEXP has_na_sexp);
extern SEXP zd_gower_dist_cross_R(SEXP x_num, SEXP x_cat, SEXP y_num, SEXP y_cat, SEXP ranges, SEXP w_num, SEXP w_cat, SEXP num_threads_sexp, SEXP has_na_sexp);

// Define the R_CallMethodDef array
static const R_CallMethodDef CallEntries[] = {
    {"zd_euclidean_dist_self_R", (DL_FUNC) &zd_euclidean_dist_self_R, 3},
    {"zd_euclidean_dist_cross_R", (DL_FUNC) &zd_euclidean_dist_cross_R, 4},
    {"zd_gower_dist_self_R", (DL_FUNC) &zd_gower_dist_self_R, 7},
    {"zd_gower_dist_cross_R", (DL_FUNC) &zd_gower_dist_cross_R, 9},
    {NULL, NULL, 0}
};

void R_init_zigdist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
