import ctypes, os, sys
from ctypes import c_double, c_int, c_size_t, POINTER, byref

# Resolve artifacts path relative to this file
_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_art = os.path.join(_root, "artifacts")

if sys.platform.startswith("win"):
    libpath = os.path.join(_art, "windows-x86_64", "msut.dll")
elif sys.platform == "darwin":
    libpath = os.path.join(_art, f"macos-{os.uname().machine}", "libmsut.dylib")
else:
    libpath = os.path.join(_art, f"linux-{os.uname().machine}", "libmsut.so")

_lib = ctypes.CDLL(libpath)

# Signatures
_lib.add_i32.argtypes = (c_int, c_int)
_lib.add_i32.restype  = c_int

_lib.sum_f64.argtypes = (POINTER(c_double), c_size_t, POINTER(c_double))
_lib.sum_f64.restype  = c_int

def add(a: int, b: int) -> int:
    return _lib.add_i32(int(a), int(b))

def sum_array(xs) -> float:
    xs = list(xs)
    arr = (c_double * len(xs))(*map(float, xs))
    out = c_double()
    rc = _lib.sum_f64(arr, len(xs), byref(out))
    if rc != 0:
        raise RuntimeError(f"sum_f64 failed with code {rc}")
    return float(out.value)
