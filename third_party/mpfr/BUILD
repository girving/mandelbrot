# GMP

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

genrule(
    name = "header",
    outs = ["mpfr.h"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/include/mpfr.h $(OUTS)",
        "//conditions:default": "\n".join([
            "if [ -e /usr/include/mpfr.h ]; then",
            "  ln -s /usr/include/mpfr.h $(OUTS)",
            "elif [ -e /usr/include/x86_64-linux-gnu/mpfr.h ]; then",
            "  ln -s /usr/include/x86_64-linux-gnu/mpfr.h $(OUTS)",
            "fi",
        ]),
    }),
)

genrule(
    name = "dylib",
    outs = ["libmpfr-link.dylib"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/lib/libmpfr.dylib $(OUTS)",
        "//conditions:default": "touch $(OUTS)",
    }),
)

genrule(
    name = "so",
    outs = ["libmpfr-link.so"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "touch $(OUTS)",
        "//conditions:default": "\n".join([
            "if [ -e /usr/lib/libmpfr.so ]; then",
            "  ln -s /usr/lib/libmpfr.so $(OUTS)",
            "elif [ -e /usr/lib/x86_64-linux-gnu/libmpfr.so ]; then",
            "  ln -s /usr/lib/x86_64-linux-gnu/libmpfr.so $(OUTS)",
            "fi",
        ]),
    }),
)

cc_library(
    name = "mpfr",
    srcs = select({
        "@bazel_tools//src/conditions:darwin": ["libmpfr-link.dylib"],
        "//conditions:default": ["libmpfr-link.so"],
    }),
    hdrs = ["mpfr.h"],
    includes = ["."],
)
