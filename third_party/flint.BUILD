# Flint

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

genrule(
    name = "fmpz-conversions",
    srcs = ["fmpz-conversions-single.in"],
    outs = ["fmpz-conversions.h"],
    cmd = "cp $(SRCS) $(OUTS)",
)

genrule(
    name = "fft_tuning",
    srcs = ["fft_tuning64.in"],
    outs = ["fft_tuning.h"],
    cmd = "cp $(SRCS) $(OUTS)",
)

genrule(
    name = "CPimport",
    srcs = ["qadic/CPimport.txt"],
    outs = ["CPimport.h"],
    cmd = "sed \"s/ /,/g;s/.*/&,/g\" $(SRCS) > $(OUTS)",
)

cc_library(
    name = "flint",
    srcs = glob([
        "*.c",
        "*/*.c",
    ], exclude = [
        "examples/*.c",
        "profile/*.c",
        "test/*.c",
        "gettimeofday.c",  # Windows only
    ]) + [
        "CPimport.h",
        "fft_tuning.h",
        "fmpz/link/fmpz_single.c",
    ],
    hdrs = glob([
        "*.h",
        "*_templates/*.c",
        "fmpz_lll/*.c",
    ]) + [
        "fmpz-conversions.h",
    ],
    include_prefix = "flint",
    copts = [
        "-Wall",
        "-Werror",
        "-Wno-shift-negative-value",
    ] + select({
        "@bazel_tools//src/conditions:darwin": [],
        "//conditions:default": ["-Wno-unused-but-set-variable"],
    }),
    deps = [
        "@mpfr//:mpfr",
        "@mandelbrot//third_party:flint-config",
    ],
)
