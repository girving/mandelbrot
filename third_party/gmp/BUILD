# GMP

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

genrule(
    name = "header",
    outs = ["gmp.h"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/include/gmp.h $(OUTS)",
        "//conditions:default": "\n".join([
            "if [ -e /usr/include/gmp.h ]; then",
            "  ln -s /usr/include/gmp.h $(OUTS)",
            "elif [ -e /usr/include/x86_64-linux-gnu/gmp.h ]; then",
            "  ln -s /usr/include/x86_64-linux-gnu/gmp.h $(OUTS)",
            "fi",
        ]),
    }),
)

genrule(
    name = "dylib",
    outs = ["libgmp-link.dylib"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/lib/libgmp.dylib $(OUTS)",
        "//conditions:default": "touch $(OUTS)",
    }),
)

genrule(
    name = "so",
    outs = ["libgmp-link.so"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "touch $(OUTS)",
        "//conditions:default": "\n".join([
            "if [ -e /usr/lib/libgmp.so ]; then",
            "  ln -s /usr/lib/libgmp.so $(OUTS)",
            "elif [ -e /usr/lib/x86_64-linux-gnu/libgmp.so ]; then",
            "  ln -s /usr/lib/x86_64-linux-gnu/libgmp.so $(OUTS)",
            "fi",
        ]),
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
