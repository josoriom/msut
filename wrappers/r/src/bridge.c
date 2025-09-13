#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#if defined(_WIN32)
#include <windows.h>
#define DLIB HMODULE
#define DLOPEN(p) LoadLibraryA(p)
#define DLSYM(h, s) GetProcAddress(h, s)
#define DLCLOSE(h) FreeLibrary(h)
static const char *last_err = "LoadLibrary/GetProcAddress failed";
#else
#include <dlfcn.h>
#define DLIB void *
#define DLOPEN(p) dlopen(p, RTLD_NOW | RTLD_GLOBAL)
#define DLSYM(h, s) dlsym(h, s)
#define DLCLOSE(h) dlclose(h)
static const char *last_err = NULL;
#endif

typedef struct
{
  unsigned char *ptr;
  size_t len;
} Buf;

typedef struct
{
  double integral_threshold;
  double intensity_threshold;
  int32_t width_threshold;
  double noise;
  int32_t auto_noise;
  int32_t allow_overlap;
  int32_t window_size;
  int32_t sn_ratio;
} CPeakPOptions;

typedef int32_t (*fn_parse_mzml)(const unsigned char *, size_t, Buf *);
typedef int32_t (*fn_bin_to_json)(const unsigned char *, size_t, Buf *);
typedef int32_t (*fn_get_peak)(const double *, const float *, size_t, double, double, const CPeakPOptions *, Buf *);
typedef int32_t (*fn_calculate_eic)(const unsigned char *, size_t, const unsigned char *, size_t, double, double, double, double, Buf *, Buf *);
typedef float (*fn_find_noise_level)(const float *, size_t);
typedef int32_t (*fn_get_peaks_from_eic)(const unsigned char *, size_t, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, double, double, const CPeakPOptions *, size_t, Buf *);
typedef int32_t (*fn_get_peaks_from_chrom)(const unsigned char *, size_t, const uint32_t *, const double *, const double *, size_t, const CPeakPOptions *, size_t, Buf *);
typedef int32_t (*fn_find_peaks)(const double *, const float *, size_t, const CPeakPOptions *, Buf *);
typedef void (*fn_free_)(unsigned char *, size_t);

typedef struct
{
  fn_parse_mzml parse_mzml;
  fn_bin_to_json bin_to_json;
  fn_get_peak get_peak;
  fn_calculate_eic calculate_eic;
  fn_find_noise_level find_noise_level;
  fn_get_peaks_from_eic C_get_peaks_from_eic;
  fn_get_peaks_from_chrom C_get_peaks_from_chrom;
  fn_find_peaks find_peaks;
  fn_free_ free_;
} abi_type;

static DLIB abi_handle = NULL;
abi_type ABI = (abi_type){0};

static int resolve_required(void **fn, const char *name)
{
  *fn = DLSYM(abi_handle, name);
  return *fn ? 0 : -1;
}
static void resolve_optional2(void **fn, const char *s1, const char *s2)
{
  *fn = DLSYM(abi_handle, s1);
  if (!*fn && s2)
    *fn = DLSYM(abi_handle, s2);
}

int abi_load(const char *path, const char **err)
{
  if (abi_handle)
  {
    DLCLOSE(abi_handle);
    abi_handle = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));

  abi_handle = DLOPEN(path);
#if !defined(_WIN32)
  if (!abi_handle)
    last_err = dlerror();
#endif
  if (!abi_handle)
  {
    if (err)
      *err = last_err;
    return -1;
  }

  if (resolve_required((void **)&ABI.parse_mzml, "parse_mzml"))
    goto fail;
  if (resolve_required((void **)&ABI.bin_to_json, "bin_to_json"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peak, "get_peak"))
    goto fail;
  if (resolve_required((void **)&ABI.calculate_eic, "calculate_eic"))
    goto fail;

  resolve_optional2((void **)&ABI.find_noise_level, "find_noise_level", NULL);
  resolve_optional2((void **)&ABI.C_get_peaks_from_eic, "C_get_peaks_from_eic", "get_peaks_from_eic");
  resolve_optional2((void **)&ABI.C_get_peaks_from_chrom, "C_get_peaks_from_chrom", "get_peaks_from_chrom");
  resolve_optional2((void **)&ABI.find_peaks, "find_peaks", "C_find_peaks");

  ABI.free_ = (fn_free_)DLSYM(abi_handle, "free_");
  if (!ABI.free_)
    goto fail;
  return 0;

fail:
  if (abi_handle)
  {
    DLCLOSE(abi_handle);
    abi_handle = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
  if (err)
    *err = last_err;
  return -1;
}

void abi_unload(void)
{
  if (abi_handle)
  {
    DLCLOSE(abi_handle);
    abi_handle = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
}

static void die_code(const char *fname, int code)
{
  const char *msg = "unknown error";
  if (code == 0)
    return;
  if (code == 1)
    msg = "invalid arguments";
  else if (code == 2)
    msg = "panic inside Rust";
  else if (code == 4)
    msg = "parse error";
  error("msut/%s failed: %s (code=%d)", fname, msg, code);
}

static SEXP mk_string_len(const unsigned char *ptr, size_t len)
{
  SEXP s = PROTECT(Rf_ScalarString(Rf_mkCharLenCE((const char *)ptr, (int)len, CE_UTF8)));
  UNPROTECT(1);
  return s;
}

#define REQUIRE_BOUND(ptr, name)                                         \
  do                                                                     \
  {                                                                      \
    if ((ptr) == NULL)                                                   \
      error("msut: symbol %s is not bound; did .onLoad() run?", (name)); \
  } while (0)

static SEXP list_get(SEXP lst, const char *name)
{
  if (TYPEOF(lst) != VECSXP)
    return R_NilValue;
  SEXP names = Rf_getAttrib(lst, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP)
    return R_NilValue;
  R_xlen_t n = XLENGTH(lst);
  for (R_xlen_t i = 0; i < n; i++)
  {
    SEXP nm = STRING_ELT(names, i);
    if (nm == R_NilValue)
      continue;
    if (strcmp(CHAR(nm), name) == 0)
      return VECTOR_ELT(lst, i);
  }
  return R_NilValue;
}

static int fill_options(SEXP opts, CPeakPOptions *out)
{
  if (opts == R_NilValue || TYPEOF(opts) != VECSXP || XLENGTH(opts) == 0)
    return 0;
  out->integral_threshold = NAN;
  out->intensity_threshold = NAN;
  out->width_threshold = 0;
  out->noise = NAN;
  out->auto_noise = 0;
  out->allow_overlap = 0;
  out->window_size = 0;
  out->sn_ratio = 0;

  SEXP v = R_NilValue;
  v = list_get(opts, "integral_threshold");
  if (v != R_NilValue)
    out->integral_threshold = asReal(v);
  v = list_get(opts, "intensity_threshold");
  if (v != R_NilValue)
    out->intensity_threshold = asReal(v);
  v = list_get(opts, "width_threshold");
  if (v != R_NilValue)
    out->width_threshold = (int32_t)asInteger(v);
  v = list_get(opts, "noise");
  if (v != R_NilValue)
    out->noise = asReal(v);
  v = list_get(opts, "auto_noise");
  if (v != R_NilValue)
    out->auto_noise = (int32_t)asLogical(v);
  v = list_get(opts, "allow_overlap");
  if (v != R_NilValue)
    out->allow_overlap = (int32_t)asLogical(v);
  v = list_get(opts, "window_size");
  if (v != R_NilValue)
    out->window_size = (int32_t)asInteger(v);
  v = list_get(opts, "sn_ratio");
  if (v != R_NilValue)
    out->sn_ratio = (int32_t)asInteger(v);
  return 1;
}

SEXP C_bind_rust(SEXP path_)
{
  if (TYPEOF(path_) != STRSXP || LENGTH(path_) != 1)
    error("path");
  const char *path = CHAR(STRING_ELT(path_, 0));
  const char *err = NULL;
  if (abi_load(path, &err) != 0)
    error("dlopen failed: %s", err ? err : "unknown");
  return R_NilValue;
}

SEXP C_parse_mzml(SEXP data)
{
  if (TYPEOF(data) != RAWSXP)
    error("data");
  REQUIRE_BOUND(ABI.parse_mzml, "parse_mzml");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.parse_mzml((const unsigned char *)RAW(data), (size_t)XLENGTH(data), &out);
  die_code("parse_mzml", code);
  SEXP res = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)out.len));
  memcpy(RAW(res), out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  UNPROTECT(1);
  return res;
}

SEXP C_bin_to_json(SEXP bin)
{
  if (TYPEOF(bin) != RAWSXP)
    error("bin");
  REQUIRE_BOUND(ABI.bin_to_json, "bin_to_json");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.bin_to_json((const unsigned char *)RAW(bin), (size_t)XLENGTH(bin), &out);
  die_code("bin_to_json", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_get_peak(SEXP x, SEXP y, SEXP rt, SEXP range, SEXP options)
{
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP)
    error("numeric");
  if (XLENGTH(x) != XLENGTH(y) || XLENGTH(x) < 3)
    error("length");
  REQUIRE_BOUND(ABI.get_peak, "get_peak");
  REQUIRE_BOUND(ABI.free_, "free_");

  R_xlen_t n = XLENGTH(y);
  float *fy = (float *)R_alloc((size_t)n, sizeof(float));
  for (R_xlen_t i = 0; i < n; i++)
    fy[i] = (float)REAL(y)[i];

  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  if (fill_options(options, &opts))
    opt_ptr = &opts;

  Buf out = (Buf){0};
  int code = ABI.get_peak(REAL(x), fy, (size_t)n, asReal(rt), asReal(range), opt_ptr, &out);
  die_code("get_peak", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_get_peaks_from_eic(SEXP bin, SEXP rts, SEXP mzs, SEXP ranges, SEXP ids, SEXP from_left, SEXP to_right, SEXP options, SEXP cores)
{
  if (TYPEOF(bin) != RAWSXP || TYPEOF(rts) != REALSXP || TYPEOF(mzs) != REALSXP || TYPEOF(ranges) != REALSXP)
    error("bad args");
  if (!(XLENGTH(rts) == XLENGTH(mzs) && XLENGTH(mzs) == XLENGTH(ranges)))
    error("length mismatch");
  REQUIRE_BOUND(ABI.C_get_peaks_from_eic, "get_peaks_from_eic");
  REQUIRE_BOUND(ABI.free_, "free_");

  R_xlen_t n = XLENGTH(rts);
  uint32_t *offs = NULL, *lens = NULL;
  unsigned char *ids_buf = NULL;
  size_t ids_len = 0;

  if (ids != R_NilValue)
  {
    if (TYPEOF(ids) != STRSXP)
      error("ids must be character");
    offs = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
    lens = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
    size_t total = 0;
    for (R_xlen_t i = 0; i < n; i++)
    {
      SEXP s = STRING_ELT(ids, i);
      if (s != R_NilValue)
        total += (size_t)LENGTH(s);
    }
    ids_buf = (unsigned char *)R_alloc(total, 1);
    ids_len = total;
    size_t cur = 0;
    for (R_xlen_t i = 0; i < n; i++)
    {
      SEXP s = STRING_ELT(ids, i);
      if (s == R_NilValue)
      {
        offs[i] = 0;
        lens[i] = 0;
      }
      else
      {
        size_t L = (size_t)LENGTH(s);
        offs[i] = (uint32_t)cur;
        lens[i] = (uint32_t)L;
        memcpy(ids_buf + cur, (const unsigned char *)CHAR(s), L);
        cur += L;
      }
    }
  }

  size_t ncores = (cores == R_NilValue) ? 1 : (size_t)asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  if (fill_options(options, &opts))
    opt_ptr = &opts;

  Buf out = (Buf){0};
  int code = ABI.C_get_peaks_from_eic(
      (const unsigned char *)RAW(bin), (size_t)XLENGTH(bin),
      REAL(rts), REAL(mzs), REAL(ranges),
      (const uint32_t *)offs, (const uint32_t *)lens,
      (const unsigned char *)ids_buf, (size_t)ids_len,
      (size_t)n, asReal(from_left), asReal(to_right),
      opt_ptr, ncores, &out);
  die_code("get_peaks_from_eic", code);

  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_get_peaks_from_chrom(SEXP bin, SEXP idxs, SEXP rts, SEXP ranges, SEXP options, SEXP cores)
{
  if (TYPEOF(bin) != RAWSXP)
    error("bin");
  if (TYPEOF(rts) != REALSXP)
    error("rt");
  if (TYPEOF(ranges) != REALSXP)
    error("range");
  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(ranges) != n || XLENGTH(idxs) != n)
    error("length");
  REQUIRE_BOUND(ABI.C_get_peaks_from_chrom, "get_peaks_from_chrom");
  REQUIRE_BOUND(ABI.free_, "free_");

  uint32_t *uidx = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  if (TYPEOF(idxs) == INTSXP)
  {
    int *ix = INTEGER(idxs);
    for (R_xlen_t i = 0; i < n; i++)
    {
      int v = ix[i];
      uidx[i] = (v == NA_INTEGER || v < 0) ? UINT32_MAX : (uint32_t)v;
    }
  }
  else if (TYPEOF(idxs) == REALSXP)
  {
    double *dx = REAL(idxs);
    for (R_xlen_t i = 0; i < n; i++)
    {
      double v = dx[i];
      uidx[i] = (!R_finite(v) || v < 0) ? UINT32_MAX : (uint32_t)v;
    }
  }
  else
    error("idx must be integer/numeric");

  size_t ncores = (cores == R_NilValue) ? 1 : (size_t)asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  if (fill_options(options, &opts))
    opt_ptr = &opts;

  Buf out = (Buf){0};
  int code = ABI.C_get_peaks_from_chrom(
      (const unsigned char *)RAW(bin), (size_t)XLENGTH(bin),
      uidx, REAL(rts), REAL(ranges), (size_t)n,
      opt_ptr, ncores, &out);
  die_code("get_peaks_from_chrom", code);

  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_calculate_eic(SEXP bin, SEXP targets, SEXP from, SEXP to, SEXP ppm_tol, SEXP mz_tol)
{
  if (TYPEOF(bin) != RAWSXP)
    error("bin");
  if (TYPEOF(targets) != STRSXP || LENGTH(targets) != 1)
    error("targets");
  REQUIRE_BOUND(ABI.calculate_eic, "calculate_eic");
  REQUIRE_BOUND(ABI.free_, "free_");

  const char *t = CHAR(STRING_ELT(targets, 0));
  size_t tlen = strlen(t);

  Buf bx = (Buf){0}, by = (Buf){0};
  int code = ABI.calculate_eic(
      (const unsigned char *)RAW(bin), (size_t)XLENGTH(bin),
      (const unsigned char *)t, tlen,
      asReal(from), asReal(to),
      asReal(ppm_tol), asReal(mz_tol),
      &bx, &by);
  die_code("calculate_eic", code);

  size_t nx = bx.len / 8;
  SEXP Rx = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nx));
  memcpy(REAL(Rx), bx.ptr, bx.len);

  size_t ny = by.len / 4;
  SEXP Ry = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)ny));
  for (size_t i = 0; i < ny; i++)
  {
    float fv;
    memcpy(&fv, by.ptr + 4 * i, 4);
    REAL(Ry)
    [i] = (double)fv;
  }

  ABI.free_(bx.ptr, bx.len);
  ABI.free_(by.ptr, by.len);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, Rx);
  SET_VECTOR_ELT(out, 1, Ry);
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("x"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("y"));
  Rf_setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(4);
  return out;
}

SEXP C_find_peaks_json(SEXP x, SEXP y, SEXP options)
{
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP)
    error("numeric");
  if (XLENGTH(x) != XLENGTH(y) || XLENGTH(x) < 3)
    error("length");
  REQUIRE_BOUND(ABI.find_peaks, "find_peaks");
  REQUIRE_BOUND(ABI.free_, "free_");

  R_xlen_t n = XLENGTH(y);
  float *fy = (float *)R_alloc((size_t)n, sizeof(float));
  for (R_xlen_t i = 0; i < n; i++)
    fy[i] = (float)REAL(y)[i];

  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  if (fill_options(options, &opts))
    opt_ptr = &opts;

  Buf out = (Buf){0};
  int code = ABI.find_peaks(REAL(x), fy, (size_t)n, opt_ptr, &out);
  die_code("find_peaks", code);

  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}
