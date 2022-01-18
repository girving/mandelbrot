package(default_visibility = ["//visibility:public"])
load("//:mandelbrot.bzl", "cc_tests")

cc_library(
    name = "known",
    hdrs = ["known.h"],
    copts = ["-Wall", "-Werror"],
)

cc_library(
    name = "mandelbrot",
    srcs = [
      "arb-cc.h",
      "arb-cc.cc",
      "arf-cc.h",
      "arf-cc.cc",
      "poly.h",
      "poly.cc",
      "print.h",
      "rand.h",
      "wall_time.h",
    ],
    copts = ["-std=c++2a", "-Wall", "-Werror"],
    deps = [
        ":known",
        "@arb//:arb",
        "@tinyformat//:tinyformat",
    ],
)

cc_binary(
    name = "area",
    srcs = ["area.cc"],
    copts = ["-std=c++2a", "-Wall", "-Werror"],
    deps = [
        ":known",
        ":mandelbrot",
    ],
)


cc_tests(
    names = ["poly_test"],
    deps = [":mandelbrot"],
)
