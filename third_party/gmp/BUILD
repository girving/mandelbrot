# GMP

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

genrule(
    name = "header",
    outs = ["gmp.h"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/include/gmp.h $(OUTS)",
        "//conditions:default": "ln -s /usr/include/gmp.h $(OUTS)",
    }),
)

genrule(
    name = "dylib",
    outs = ["libgmp-link.dylib"],
    cmd = "ln -s `brew --prefix`/lib/libgmp.dylib $(OUTS)",
)

genrule(
    name = "so",
    outs = ["libgmp-link.so"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "touch $(OUTS)",
        "//conditions:default": "ln -s `brew --prefix`/lib/libgmp.so $(OUTS)",
    }),
)

cc_library(
    name = "gmp",
    srcs = select({
        "@bazel_tools//src/conditions:darwin": ["libgmp-link.dylib"],
        "//conditions:default": ["libgmp-link.so"],
    }),
    hdrs = ["gmp.h"],
    includes = ["."],
)
