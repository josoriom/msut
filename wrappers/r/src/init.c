#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP C_bind_rust(SEXP path);
SEXP C_parse_mzml(SEXP data);
SEXP C_bin_to_json(SEXP bin);
SEXP C_get_peak(SEXP x, SEXP y, SEXP rt, SEXP range, SEXP options);
SEXP C_get_peaks_from_eic(SEXP bin, SEXP rts, SEXP mzs, SEXP ranges, SEXP ids, SEXP from_left, SEXP to_right, SEXP options, SEXP cores);
SEXP C_get_peaks_from_chrom(SEXP bin, SEXP idxs, SEXP rts, SEXP ranges, SEXP options, SEXP cores);
SEXP C_calculate_eic(SEXP bin, SEXP targets, SEXP from, SEXP to, SEXP ppm_tol, SEXP mz_tol);
SEXP C_find_peaks_json(SEXP x, SEXP y, SEXP options);

static const R_CallMethodDef CallEntries[] = {
    {"C_bind_rust", (DL_FUNC)&C_bind_rust, 1},
    {"C_parse_mzml", (DL_FUNC)&C_parse_mzml, 1},
    {"C_bin_to_json", (DL_FUNC)&C_bin_to_json, 1},
    {"C_get_peak", (DL_FUNC)&C_get_peak, 5},
    {"C_get_peaks_from_eic", (DL_FUNC)&C_get_peaks_from_eic, 9},
    {"C_get_peaks_from_chrom", (DL_FUNC)&C_get_peaks_from_chrom, 6},
    {"C_calculate_eic", (DL_FUNC)&C_calculate_eic, 6},
    {"C_find_peaks_json", (DL_FUNC)&C_find_peaks_json, 3},
    {NULL, NULL, 0}};

void R_init_msut(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
