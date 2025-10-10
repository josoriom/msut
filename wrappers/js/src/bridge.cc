#if __has_include(<napi.h>)
#include <napi.h>
#else
#include "node-addon-api/napi.h"
#endif

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <string>
#include <vector>

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
  int32_t auto_baseline;
  int32_t baseline_window;
  int32_t baseline_window_factor;
  int32_t allow_overlap;
  int32_t window_size;
  double sn_ratio;
} CPeakPOptions;

static_assert(sizeof(CPeakPOptions) == 64, "CPeakPOptions must be 64 bytes");

typedef int32_t (*fn_parse_mzml)(const unsigned char *, size_t, Buf *);
typedef int32_t (*fn_bin_to_json)(const unsigned char *, size_t, Buf *);
typedef int32_t (*fn_get_peak)(
    const double *, const double *, size_t, double, double, const CPeakPOptions *, Buf *);
typedef int32_t (*fn_calculate_eic)(
    const unsigned char *, size_t,
    double,
    double, double, double, double,
    Buf *, Buf *);
typedef double (*fn_find_noise_level)(const double *, size_t);
typedef int32_t (*fn_get_peaks_from_eic)(
    const unsigned char *, size_t,
    const double *, const double *, const double *,
    const uint32_t *, const uint32_t *, const unsigned char *, size_t,
    size_t, double, double, const CPeakPOptions *, size_t, Buf *);
typedef int32_t (*fn_get_peaks_from_chrom)(
    const unsigned char *, size_t,
    const uint32_t *, const double *, const double *, size_t,
    const CPeakPOptions *, size_t, Buf *);
typedef int32_t (*fn_find_peaks)(
    const double *, const double *, size_t, const CPeakPOptions *, Buf *);
typedef int32_t (*fn_calculate_baseline)(
    const double *, size_t, int32_t, int32_t, Buf *);
typedef int32_t (*fn_find_features)(
    const unsigned char *, size_t,
    double, double,
    double, double,
    double, double, double,
    const CPeakPOptions *, int32_t,
    Buf *);
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
  fn_calculate_baseline calculate_baseline;
  fn_find_features find_features;
  fn_free_ free_;
} msabi_t;

static msabi_t ABI{};
static DLIB LIB_HANDLE = NULL;

static int resolve_required(void **fn, const char *name)
{
#if !defined(_WIN32)
  dlerror();
#endif
  *fn = DLSYM(LIB_HANDLE, name);
  if (!*fn)
  {
#if !defined(_WIN32)
    last_err = dlerror();
#endif
    return -1;
  }
  return 0;
}

static int abi_load(const char *path, const char **err)
{
  if (LIB_HANDLE)
  {
    DLCLOSE(LIB_HANDLE);
    LIB_HANDLE = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
#if !defined(_WIN32)
  dlerror();
#endif
  LIB_HANDLE = DLOPEN(path);
#if !defined(_WIN32)
  if (!LIB_HANDLE)
    last_err = dlerror();
#endif
  if (!LIB_HANDLE)
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
  if (resolve_required((void **)&ABI.C_get_peaks_from_eic, "get_peaks_from_eic"))
    goto fail;
  if (resolve_required((void **)&ABI.C_get_peaks_from_chrom, "get_peaks_from_chrom"))
    goto fail;
  if (resolve_required((void **)&ABI.find_peaks, "find_peaks"))
    goto fail;
  ABI.calculate_baseline = (fn_calculate_baseline)DLSYM(LIB_HANDLE, "calculate_baseline");
  if (!ABI.calculate_baseline)
    ABI.calculate_baseline = (fn_calculate_baseline)DLSYM(LIB_HANDLE, "calculate_baseline_v2");
  if (resolve_required((void **)&ABI.find_features, "find_features"))
    goto fail;

  ABI.find_noise_level = (fn_find_noise_level)DLSYM(LIB_HANDLE, "find_noise_level");
  ABI.free_ = (fn_free_)DLSYM(LIB_HANDLE, "free_");
  if (!ABI.free_)
    goto fail;

  return 0;

fail:
  if (LIB_HANDLE)
  {
    DLCLOSE(LIB_HANDLE);
    LIB_HANDLE = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
  if (err)
    *err = last_err;
  return -1;
}

static void abi_unload(void)
{
  if (LIB_HANDLE)
  {
    DLCLOSE(LIB_HANDLE);
    LIB_HANDLE = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
}

static const char *CodeMessage(int32_t code)
{
  if (code == 0)
    return "ok";
  if (code == 1)
    return "invalid arguments";
  if (code == 2)
    return "panic inside Rust";
  if (code == 4)
    return "parse error";
  return "unknown";
}

static Napi::Buffer<uint8_t> TakeBuffer(Napi::Env env, Buf *buf)
{
  Napi::Buffer<uint8_t> out = Napi::Buffer<uint8_t>::Copy(env, buf->ptr, buf->len);
  if (ABI.free_)
    ABI.free_(buf->ptr, buf->len);
  buf->ptr = nullptr;
  buf->len = 0;
  return out;
}

static const CPeakPOptions *ReadOptionsBuf(Napi::Value value, CPeakPOptions *out)
{
  if (value.IsUndefined() || value.IsNull())
    return nullptr;
  if (!value.IsBuffer())
    return nullptr;
  Napi::Buffer<uint8_t> buf = value.As<Napi::Buffer<uint8_t>>();
  if (buf.Length() != sizeof(CPeakPOptions))
    return nullptr;
  memcpy(out, buf.Data(), sizeof(CPeakPOptions));
  return out;
}

static void ThrowIfMissing(Napi::Env env, void *fn_ptr, const char *name)
{
  if (fn_ptr == nullptr)
  {
    std::string msg = "native symbol not exported: ";
    msg += name;
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
  }
}

static Napi::Value Bind(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 1 || !info[0].IsString())
  {
    Napi::TypeError::New(env, "expected: path string").ThrowAsJavaScriptException();
    return env.Undefined();
  }
  std::string path = info[0].As<Napi::String>();
  const char *err = nullptr;
  if (abi_load(path.c_str(), &err) != 0)
  {
    std::string msg = "dlopen failed: ";
    if (err)
      msg += err;
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  return env.Undefined();
}

static Napi::Value ParseMzML(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.parse_mzml, "parse_mzml");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");
  Napi::Buffer<uint8_t> input = info[0].As<Napi::Buffer<uint8_t>>();
  Buf out = {nullptr, 0};
  int32_t rc = ABI.parse_mzml(input.Data(), (size_t)input.Length(), &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "parse_mzml: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  return TakeBuffer(env, &out);
}

static Napi::Value BinToJson(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.bin_to_json, "bin_to_json");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");
  Napi::Buffer<uint8_t> bin = info[0].As<Napi::Buffer<uint8_t>>();
  Buf out = {nullptr, 0};
  int32_t rc = ABI.bin_to_json(bin.Data(), (size_t)bin.Length(), &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "bin_to_json: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  std::string json_text((const char *)out.ptr, out.len);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  return Napi::String::New(env, json_text);
}

static Napi::Value GetPeak(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.get_peak, "get_peak");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  Napi::Float64Array x_arr = info[0].As<Napi::Float64Array>();
  Napi::Float64Array y_arr = info[1].As<Napi::Float64Array>();
  double target_rt = info[2].As<Napi::Number>().DoubleValue();
  double rt_range = info[3].As<Napi::Number>().DoubleValue();

  CPeakPOptions opts;
  const CPeakPOptions *p_opts = nullptr;
  if (info.Length() > 4)
    p_opts = ReadOptionsBuf(info[4], &opts);

  double *x_ptr = (double *)((uint8_t *)x_arr.ArrayBuffer().Data() + x_arr.ByteOffset());
  double *y_ptr = (double *)((uint8_t *)y_arr.ArrayBuffer().Data() + y_arr.ByteOffset());
  size_t n = x_arr.ElementLength();

  Buf out = {nullptr, 0};
  int32_t rc = ABI.get_peak(x_ptr, y_ptr, n, target_rt, rt_range, p_opts, &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "get_peak: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  std::string json_text((const char *)out.ptr, out.len);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  return Napi::String::New(env, json_text);
}

static Napi::Value CalculateEic(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.calculate_eic, "calculate_eic");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  Napi::Buffer<uint8_t> bin = info[0].As<Napi::Buffer<uint8_t>>();
  double targets = info[1].As<Napi::Number>().DoubleValue();
  double from_rt = info[2].As<Napi::Number>().DoubleValue();
  double to_rt = info[3].As<Napi::Number>().DoubleValue();
  double ppm_tol = info[4].As<Napi::Number>().DoubleValue();
  double mz_tol = info[5].As<Napi::Number>().DoubleValue();

  Buf x_buf = {nullptr, 0};
  Buf y_buf = {nullptr, 0};
  int32_t rc = ABI.calculate_eic(
      bin.Data(), (size_t)bin.Length(),
      targets,
      from_rt, to_rt, ppm_tol, mz_tol,
      &x_buf, &y_buf);
  if (rc != 0)
  {
    if (x_buf.ptr && ABI.free_)
      ABI.free_(x_buf.ptr, x_buf.len);
    if (y_buf.ptr && ABI.free_)
      ABI.free_(y_buf.ptr, y_buf.len);
    std::string msg = "calculate_eic: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }

  size_t nx = x_buf.len / 8, ny = y_buf.len / 8;
  Napi::ArrayBuffer abx = Napi::ArrayBuffer::New(env, nx * 8);
  Napi::ArrayBuffer aby = Napi::ArrayBuffer::New(env, ny * 8);
  memcpy(abx.Data(), x_buf.ptr, nx * 8);
  memcpy(aby.Data(), y_buf.ptr, ny * 8);
  if (ABI.free_)
  {
    ABI.free_(x_buf.ptr, x_buf.len);
    ABI.free_(y_buf.ptr, y_buf.len);
  }
  Napi::Float64Array X = Napi::Float64Array::New(env, nx, abx, 0);
  Napi::Float64Array Y = Napi::Float64Array::New(env, ny, aby, 0);
  Napi::Object out = Napi::Object::New(env);
  out.Set("x", X);
  out.Set("y", Y);
  return out;
}

static Napi::Value FindNoiseLevel(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.find_noise_level, "find_noise_level");
  Napi::Float64Array y_arr = info[0].As<Napi::Float64Array>();
  double *y_ptr = (double *)((uint8_t *)y_arr.ArrayBuffer().Data() + y_arr.ByteOffset());
  double value = ABI.find_noise_level(y_ptr, y_arr.ElementLength());
  return Napi::Number::New(env, value);
}

static Napi::Value GetPeaksFromEic(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.C_get_peaks_from_eic, "get_peaks_from_eic");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  Napi::Buffer<uint8_t> bin = info[0].As<Napi::Buffer<uint8_t>>();
  Napi::Float64Array rts_arr = info[1].As<Napi::Float64Array>();
  Napi::Float64Array mzs_arr = info[2].As<Napi::Float64Array>();
  Napi::Float64Array rng_arr = info[3].As<Napi::Float64Array>();

  const double *rts = (const double *)((uint8_t *)rts_arr.ArrayBuffer().Data() + rts_arr.ByteOffset());
  const double *mzs = (const double *)((uint8_t *)mzs_arr.ArrayBuffer().Data() + mzs_arr.ByteOffset());
  const double *rng = (const double *)((uint8_t *)rng_arr.ArrayBuffer().Data() + rng_arr.ByteOffset());
  size_t count = rts_arr.ElementLength();

  const uint32_t *offs_ptr = nullptr;
  const uint32_t *lens_ptr = nullptr;
  const unsigned char *ids_buf_ptr = nullptr;
  size_t ids_buf_len = 0;

  Napi::Buffer<uint32_t> offs_js;
  Napi::Buffer<uint32_t> lens_js;
  Napi::Buffer<unsigned char> ids_js;

  if (info.Length() > 4 && !info[4].IsUndefined() && !info[4].IsNull())
  {
    Napi::Array ids = info[4].As<Napi::Array>();
    std::vector<uint32_t> offs(count, 0);
    std::vector<uint32_t> lens(count, 0);
    size_t total = 0;
    std::vector<std::string> tmp;
    tmp.reserve(count);
    for (size_t i = 0; i < count; i++)
    {
      Napi::Value v = ids.Get((uint32_t)i);
      if (v.IsString())
      {
        std::string s = v.As<Napi::String>().Utf8Value();
        total += s.size();
        tmp.push_back(std::move(s));
      }
      else
      {
        tmp.emplace_back();
      }
    }
    std::vector<unsigned char> ids_buf;
    ids_buf.resize(total);
    size_t cur = 0;
    for (size_t i = 0; i < count; i++)
    {
      const std::string &s = tmp[i];
      offs[i] = (uint32_t)cur;
      lens[i] = (uint32_t)s.size();
      if (!s.empty())
      {
        memcpy(ids_buf.data() + cur, s.data(), s.size());
        cur += s.size();
      }
    }
    offs_js = Napi::Buffer<uint32_t>::Copy(env, offs.data(), offs.size());
    lens_js = Napi::Buffer<uint32_t>::Copy(env, lens.data(), lens.size());
    ids_js = Napi::Buffer<unsigned char>::Copy(env, ids_buf.data(), ids_buf.size());
    offs_ptr = (const uint32_t *)offs_js.Data();
    lens_ptr = (const uint32_t *)lens_js.Data();
    ids_buf_ptr = (const unsigned char *)ids_js.Data();
    ids_buf_len = ids_js.Length();
  }

  double from_left = info[5].As<Napi::Number>().DoubleValue();
  double to_right = info[6].As<Napi::Number>().DoubleValue();

  CPeakPOptions opts;
  const CPeakPOptions *p_opts = nullptr;
  if (info.Length() > 7)
    p_opts = ReadOptionsBuf(info[7], &opts);

  size_t cores = 1;
  if (info.Length() > 8 && info[8].IsNumber())
  {
    int64_t v = info[8].As<Napi::Number>().Int64Value();
    if (v > 0)
      cores = (size_t)v;
  }

  Buf out = {nullptr, 0};
  int32_t rc = ABI.C_get_peaks_from_eic(
      bin.Data(), (size_t)bin.Length(),
      rts, mzs, rng,
      offs_ptr, lens_ptr, ids_buf_ptr, ids_buf_len,
      count, from_left, to_right, p_opts, cores, &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "get_peaks_from_eic: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  std::string json_text((const char *)out.ptr, out.len);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  return Napi::String::New(env, json_text);
}

static Napi::Value GetPeaksFromChrom(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.C_get_peaks_from_chrom, "get_peaks_from_chrom");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  Napi::Buffer<uint8_t> bin = info[0].As<Napi::Buffer<uint8_t>>();
  Napi::Uint32Array idxs_arr = info[1].As<Napi::Uint32Array>();
  Napi::Float64Array rts_arr = info[2].As<Napi::Float64Array>();
  Napi::Float64Array rng_arr = info[3].As<Napi::Float64Array>();

  const uint32_t *idx = (const uint32_t *)((uint8_t *)idxs_arr.ArrayBuffer().Data() + idxs_arr.ByteOffset());
  const double *rts = (const double *)((uint8_t *)rts_arr.ArrayBuffer().Data() + rts_arr.ByteOffset());
  const double *rng = (const double *)((uint8_t *)rng_arr.ArrayBuffer().Data() + rng_arr.ByteOffset());
  size_t count = rts_arr.ElementLength();

  CPeakPOptions opts;
  const CPeakPOptions *p_opts = nullptr;
  if (info.Length() > 4)
    p_opts = ReadOptionsBuf(info[4], &opts);

  size_t cores = 1;
  if (info.Length() > 5 && info[5].IsNumber())
  {
    int64_t v = info[5].As<Napi::Number>().Int64Value();
    if (v > 0)
      cores = (size_t)v;
  }

  Buf out = {nullptr, 0};
  int32_t rc = ABI.C_get_peaks_from_chrom(
      bin.Data(), (size_t)bin.Length(),
      idx, rts, rng, count,
      p_opts, cores, &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "get_peaks_from_chrom: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  std::string json_text((const char *)out.ptr, out.len);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  return Napi::String::New(env, json_text);
}

static Napi::Value FindPeaks(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.find_peaks, "find_peaks");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  Napi::Float64Array x_arr = info[0].As<Napi::Float64Array>();
  Napi::Float64Array y_arr = info[1].As<Napi::Float64Array>();

  CPeakPOptions opts;
  const CPeakPOptions *p_opts = nullptr;
  if (info.Length() > 2)
    p_opts = ReadOptionsBuf(info[2], &opts);

  double *x_ptr = (double *)((uint8_t *)x_arr.ArrayBuffer().Data() + x_arr.ByteOffset());
  double *y_ptr = (double *)((uint8_t *)y_arr.ArrayBuffer().Data() + y_arr.ByteOffset());
  size_t n = x_arr.ElementLength();

  Buf out = {nullptr, 0};
  int32_t rc = ABI.find_peaks(x_ptr, y_ptr, n, p_opts, &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "find_peaks: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  std::string json_text((const char *)out.ptr, out.len);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  return Napi::String::New(env, json_text);
}

static Napi::Value CalculateBaseline(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.calculate_baseline, "calculate_baseline");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  if (info.Length() < 1 || !info[0].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: Float64Array").ThrowAsJavaScriptException();
    return env.Undefined();
  }
  Napi::Float64Array y_arr = info[0].As<Napi::Float64Array>();

  int32_t baseline_window = 0;
  int32_t baseline_window_factor = 0;

  if (info.Length() >= 2 && info[1].IsObject() && !info[1].IsBuffer() && !info[1].IsTypedArray())
  {
    Napi::Object o = info[1].As<Napi::Object>();
    if (o.Has("baselineWindow"))
      baseline_window = o.Get("baselineWindow").ToNumber().Int32Value();
    if (o.Has("baselineWindowFactor"))
      baseline_window_factor = o.Get("baselineWindowFactor").ToNumber().Int32Value();
  }
  else
  {
    if (info.Length() > 1 && info[1].IsNumber())
      baseline_window = info[1].As<Napi::Number>().Int32Value();
    if (info.Length() > 2 && info[2].IsNumber())
      baseline_window_factor = info[2].As<Napi::Number>().Int32Value();
  }

  const double *y_ptr = (const double *)((uint8_t *)y_arr.ArrayBuffer().Data() + y_arr.ByteOffset());
  size_t n = y_arr.ElementLength();

  Buf out = {nullptr, 0};
  int32_t rc = ABI.calculate_baseline(y_ptr, n, baseline_window, baseline_window_factor, &out);
  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "calculate_baseline: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  size_t ny = out.len / 8;
  Napi::ArrayBuffer aby = Napi::ArrayBuffer::New(env, ny * 8);
  memcpy(aby.Data(), out.ptr, ny * 8);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  Napi::Float64Array Y = Napi::Float64Array::New(env, ny, aby, 0);
  return Y;
}

static Napi::Value FindFeatures(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  ThrowIfMissing(env, (void *)ABI.find_features, "find_features");
  ThrowIfMissing(env, (void *)ABI.free_, "free_");

  if (info.Length() < 10)
  {
    Napi::TypeError::New(env,
                         "expected: (Buffer data, number from, number to, number eicPpm, number eicMz, "
                         "number gridStart, number gridEnd, number gridStepPpm, Buffer|null options, number cores)")
        .ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Buffer<uint8_t> data = info[0].As<Napi::Buffer<uint8_t>>();
  double from_time = info[1].As<Napi::Number>().DoubleValue();
  double to_time = info[2].As<Napi::Number>().DoubleValue();
  double eic_ppm = info[3].As<Napi::Number>().DoubleValue();
  double eic_mz = info[4].As<Napi::Number>().DoubleValue();
  double grid_start = info[5].As<Napi::Number>().DoubleValue();
  double grid_end = info[6].As<Napi::Number>().DoubleValue();
  double grid_step = info[7].As<Napi::Number>().DoubleValue();

  CPeakPOptions opts;
  const CPeakPOptions *p_opts = nullptr;
  if (!info[8].IsUndefined() && !info[8].IsNull())
  {
    if (!info[8].IsBuffer())
    {
      Napi::TypeError::New(env, "options must be a Buffer, null, or undefined")
          .ThrowAsJavaScriptException();
      return env.Undefined();
    }
    p_opts = ReadOptionsBuf(info[8], &opts);
    if (p_opts == nullptr)
    {
      Napi::TypeError::New(env, "options Buffer must be exactly 64 bytes")
          .ThrowAsJavaScriptException();
      return env.Undefined();
    }
  }

  if (!info[9].IsNumber())
  {
    Napi::TypeError::New(env, "cores must be a positive integer")
        .ThrowAsJavaScriptException();
    return env.Undefined();
  }
  int32_t cores = info[9].As<Napi::Number>().Int32Value();
  if (cores <= 0)
  {
    Napi::TypeError::New(env, "cores must be > 0").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Buf out = {nullptr, 0};
  int32_t rc = ABI.find_features(
      data.Data(), (size_t)data.Length(),
      from_time, to_time,
      eic_ppm, eic_mz,
      grid_start, grid_end, grid_step,
      p_opts, cores, &out);

  if (rc != 0)
  {
    if (out.ptr && ABI.free_)
      ABI.free_(out.ptr, out.len);
    std::string msg = "find_features: ";
    msg += CodeMessage(rc);
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }

  std::string json_text((const char *)out.ptr, out.len);
  if (ABI.free_)
    ABI.free_(out.ptr, out.len);
  return Napi::String::New(env, json_text);
}

static Napi::Object Init(Napi::Env env, Napi::Object exports)
{
  exports.Set("bind", Napi::Function::New(env, Bind));
  exports.Set("parseMzML", Napi::Function::New(env, ParseMzML));
  exports.Set("binToJson", Napi::Function::New(env, BinToJson));
  exports.Set("getPeak", Napi::Function::New(env, GetPeak));
  exports.Set("calculateEic", Napi::Function::New(env, CalculateEic));
  exports.Set("findNoiseLevel", Napi::Function::New(env, FindNoiseLevel));
  exports.Set("getPeaksFromEic", Napi::Function::New(env, GetPeaksFromEic));
  exports.Set("getPeaksFromChrom", Napi::Function::New(env, GetPeaksFromChrom));
  exports.Set("findPeaks", Napi::Function::New(env, FindPeaks));
  exports.Set("calculateBaseline", Napi::Function::New(env, CalculateBaseline));
  exports.Set("findFeatures", Napi::Function::New(env, FindFeatures));
  return exports;
}

NODE_API_MODULE(msut, Init)
