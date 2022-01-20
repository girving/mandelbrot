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
      "print.h",
      "wall_time.h",
    ],
    copts = copts,
    deps = [
        "@tinyformat//:tinyformat",
    ],
)

cc_library(
    name = "area-arb",
    srcs = [
      "arb-cc.h",
      "arb-cc.cc",
      "arf-cc.h",
      "arf-cc.cc",
      "area.h",
      "area.cc",
      "poly.h",
      "poly.cc",
      "rand.h",
    ],
    copts = copts,
    deps = [
        ":base",
        ":known",
        "@arb//:arb",
    ],
)

cc_library(
    name = "series",
    srcs = [
      "series.h",
    ],
    copts = copts,
    deps = [
        ":base",
    ],
)

cc_binary(
    name = "mandelbrot",
    srcs = ["mandelbrot.cc"],
    copts = copts,
    deps = [
        ":known",
        ":area-arb",
    ],
)

cc_tests(
    names = [
        "area_test",
        "poly_test",
    ],
    deps = [":area-arb"],
)

cc_tests(
    names = ["series_test"],
    deps = [
        ":series",
        ":area-arb",
    ],
)
