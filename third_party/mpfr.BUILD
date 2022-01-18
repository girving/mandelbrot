# MPFR

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

genrule(
    name = "mparam",
    srcs = ["src/mparam_h.in"],
    outs = ["src/mparam.h"],
    cmd = "cp $(SRCS) $(OUTS)",
)

cc_library(
    name = "mpfr",
    srcs = glob([
        "src/**/*.c",
    ], exclude = [
        "src/add1sp1_extracted.c",
        "src/jyn_asympt.c",
        "src/mul_1_extracted.c",
        "src/round_raw_generic.c",
        "src/sub1sp1_extracted.c",
    ]),
    hdrs = glob(["src/*.h"]) + [
        "src/mparam.h",
    ],
    includes = ["src"],
    copts = [
        "-Wall",
        "-Werror",
        "-Wno-tautological-constant-out-of-range-compare",
        # DEFS from mpfr/src/Makefile after configure
        "-DHAVE_INTTYPES_H=1",
        "-DHAVE_STDINT_H=1",
        "-DLT_OBJDIR=\".libs/\"",
        "-DHAVE_LITTLE_ENDIAN=1",
        "-DHAVE_CLOCK_GETTIME=1",
        "-DTIME_WITH_SYS_TIME=1",
        "-DHAVE_LOCALE_H=1",
        "-DHAVE_WCHAR_H=1",
        "-DHAVE_STDARG=1",
        "-DHAVE_SYS_TIME_H=1",
        "-DHAVE_STRUCT_LCONV_DECIMAL_POINT=1",
        "-DHAVE_STRUCT_LCONV_THOUSANDS_SEP=1",
        "-DHAVE_ALLOCA_H=1",
        "-DHAVE_ALLOCA=1",
        "-DHAVE_UINTPTR_T=1",
        "-DHAVE_VA_COPY=1",
        "-DHAVE_SETLOCALE=1",
        "-DHAVE_GETTIMEOFDAY=1",
        "-DHAVE_SIGNAL=1",
        "-DHAVE_SIGACTION=1",
        "-DHAVE_LONG_LONG=1",
        "-DHAVE_INTMAX_T=1",
        "-DMPFR_HAVE_INTMAX_MAX=1",
        "-DMPFR_HAVE_NORETURN=1",
        "-DMPFR_HAVE_BUILTIN_UNREACHABLE=1",
        "-DMPFR_HAVE_CONSTRUCTOR_ATTR=1",
        "-DMPFR_HAVE_FESETROUND=1",
        "-DHAVE_SUBNORM_DBL=1",
        "-DHAVE_SUBNORM_FLT=1",
        "-DHAVE_SIGNEDZ=1",
        "-DHAVE_ROUND=1",
        "-DHAVE_TRUNC=1",
        "-DHAVE_FLOOR=1",
        "-DHAVE_CEIL=1",
        "-DHAVE_NEARBYINT=1",
        "-DHAVE_DOUBLE_IEEE_LITTLE_ENDIAN=1",
        "-DHAVE_LDOUBLE_IEEE_EXT_LITTLE=1",
        "-DMPFR_USE_THREAD_SAFE=1",
        "-DMPFR_USE_C11_THREAD_SAFE=1",
        "-DMPFR_USE_STATIC_ASSERT=1",
        "-DHAVE_ATTRIBUTE_MODE=1",
        "-DPRINTF_L=1",
        "-DPRINTF_T=1",
        "-DPRINTF_GROUPFLAG=1",
        "-DHAVE___GMPN_SBPI1_DIVAPPR_Q=1",
        "-DHAVE___GMPN_INVERT_LIMB=1",
        "-DHAVE___GMPN_RSBLSH1_N=1",
        "-DMPFR_LONG_WITHIN_LIMB=1",
        "-DMPFR_INTMAX_WITHIN_LIMB=1",
        "-DHAVE_GETRUSAGE=1",
    ],
    deps = [
        "@mandelbrot//third_party/gmp",
    ],
)
