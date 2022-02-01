package(default_visibility = ["//visibility:public"])
load("//:mandelbrot.bzl", "cc_tests")
load("//:mandelbrot.bzl", "copts")

cc_library(
    name = "known",
    srcs = [
        "known.h",
        "known.cc",
    ],
    copts = copts,
)

cc_library(
    name = "base",
    srcs = [
      "arith.h",
      "bit.h",
      "debug.h",
      "debug.cc",
      "format.h",
      "is_interval.h",
      "noncopyable.h",
      "print.h",
      "span.h",
      "wall_time.h",
    ],
    copts = copts,
    deps = [
        "@tinyformat//:tinyformat",
    ],
)

cc_binary(
    name = "codelets",
    srcs = ["codelets.cc"],
    copts = copts,
    deps = [
        ":base",
    ],
)

genrule(
    name = "run-codelets",
    outs = [
        "gen/expansion.h",
    ],
    tools = [":codelets"],
    cmd = "$(location :codelets) $(OUTS)",
)

cc_library(
    name = "arb",
    srcs = [
      "acb_cc.h",
      "acb_cc.cc",
      "arb_cc.h",
      "arb_cc.cc",
      "arf_cc.h",
      "arf_cc.cc",
      "fmpq_cc.h",
      "fmpq_cc.cc",
      "mag_cc.h",
      "poly.h",
      "poly.cc",
      "rand.h",
    ],
    copts = copts,
    deps = [
        ":base",
        "@arb//:arb",
        "@flint//:flint",
    ],
)

cc_library(
    name = "arb_area",
    srcs = [
      "arb_area.h",
      "arb_area.cc",
    ],
    copts = copts,
    deps = [
        ":arb",
        ":base",
        ":known",
        "@arb//:arb",
    ],
)

cc_library(
    name = "area",
    srcs = [
        "area.h",
        "area.cc",
        "complex.h",
        "expansion.h",
        "expansion.cc",
        "fft.h",
        "fft.cc",
        "gen/expansion.h",
        "nearest.h",
        "nearest.cc",
        "series.h",
    ],
    copts = copts,
    deps = [
        ":arb",
        ":base",
        ":known",
    ],
)

cc_binary(
    name = "mandelbrot",
    srcs = ["mandelbrot.cc"],
    copts = copts,
    deps = [
        ":arb_area",
        ":area",
    ],
)

cc_tests(
    names = ["bit_test"],
    deps = [":base"],
)

cc_tests(
    names = ["expansion_test"],
    deps = [
        ":arb",
        ":area",
    ],
)

cc_tests(
    names = [
        "arb_area_test",
        "poly_test",
    ],
    deps = [":arb_area"],
)

cc_tests(
    names = [
        "area_test",
        "series_test",
    ],
    deps = [":area"],
)

cc_tests(
    names = ["fft_test"],
    deps = [":area"],
)
