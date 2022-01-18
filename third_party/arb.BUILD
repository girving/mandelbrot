# Arb

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

cc_library(
    name = "arb",
    srcs = glob([
        "*/*.c",
    ], exclude = [
        "examples/*.c",
    ]),
    hdrs = glob(["*.h"]),
    includes = ["."],
    copts = [
        "-Wall",
        "-Werror",
    ] + select({
        "@bazel_tools//src/conditions:darwin": [],
        "//conditions:default": ["-Wno-builtin-declaration-mismatch"],
    }),
    linkopts = ["-lpthread"],
    deps = [
        "@flint//:flint",
    ],
)
