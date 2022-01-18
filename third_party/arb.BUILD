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
    deps = [
        "@flint//:flint",
    ],
)
