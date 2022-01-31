package(default_visibility = ["//visibility:public"])
load("//:mandelbrot.bzl", "cc_tests")
load("//:mandelbrot.bzl", "copts")

cc_library(
    name = "known",
    hdrs = ["known.h"],
    copts = copts,
)

cc_library(
    name = "base",
    srcs = [
      "debug.h",
      "debug.cc",
      "format.h",
      "noncopyable.h",
      "print.h",
      "relu.h",
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
    copts = copts + ["-Wno-shorten-64-to-32"],
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
    name = "expansion",
    srcs = [
        "expansion.h",
        "expansion.cc",
        "gen/expansion.h",
    ],
    copts = copts,
    deps = [
        ":arb",
    ],
)

cc_library(
    name = "complex",
    srcs = ["complex.h"],
    copts = copts,
)

cc_library(
    name = "nearest",
    srcs = [
        "nearest.h",
        "nearest.cc",
    ],
    copts = copts,
    deps = [
        ":arb",
        ":base",
        ":complex",
        ":expansion",
    ],
)

cc_library(
    name = "series",
    srcs = [
      "area.h",
      "area.cc",
      "fft.h",
      "fft.cc",
      "series.h",
    ],
    copts = copts,
    deps = [
        ":arb",
        ":base",
        ":complex",
        ":known",
        ":nearest",
    ],
)

cc_binary(
    name = "mandelbrot",
    srcs = ["mandelbrot.cc"],
    copts = copts,
    deps = [
        ":arb_area",
        ":known",
        ":series",
    ],
)

cc_tests(
    names = ["expansion_test"],
    copts = copts + ["-Wno-shorten-64-to-32"],
    deps = [
        ":arb",
        ":nearest",
        ":expansion",
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
    deps = [
        ":series",
        ":arb_area",
    ],
)

cc_tests(
    names = ["fft_test"],
    copts = copts + ["-Wno-shorten-64-to-32"],
    deps = [
        ":series",
        ":arb_area",
    ],
)
