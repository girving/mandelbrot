# GMP

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

headers = [
    "acb.h",
    "acb_dft.h",
    "arb.h",
    "arb_poly.h",
    "arf.h",
    "dlog.h",
    "fmpr.h",
    "fmpz_extras.h",
    "mag.h",
]

genrule(
    name = "headers",
    outs = headers,
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "\n".join([
            "for h in $(OUTS); do",
            "  ln -s `brew --prefix`/include/$$(basename $$h) $$h",
            "done",
        ]),
        "//conditions:default": "\n".join([
            "if [ -e /usr/include/arb.h ]; then",
            "  src=/usr/include",
            "elif [ -e /usr/include/x86_64-linux-gnu/arb.h ]; then",
            "  src=/usr/include/x86_64-linux-gnu",
            "fi",
            "for h in $(OUTS); do",
            "  ln -s $$src/$$(basename $$h) $$h",
            "done",
        ]),
    }),
)

genrule(
    name = "dylib",
    outs = ["libarb-link.dylib"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/lib/libarb.dylib $(OUTS)",
        "//conditions:default": "touch $(OUTS)",
    }),
)

genrule(
    name = "so",
    outs = ["libarb-link.so"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "touch $(OUTS)",
        "//conditions:default": "\n".join([
            "if [ -e /usr/lib/libflint-arb.so ]; then",
            "  ln -s /usr/lib/libflint-arb.so $(OUTS)",
            "elif [ -e /usr/lib/x86_64-linux-gnu/libflint-arb.so ]; then",
            "  ln -s /usr/lib/x86_64-linux-gnu/libflint-arb.so $(OUTS)",
            "fi",
        ]),
    }),
)

cc_library(
    name = "arb",
    srcs = select({
        "@bazel_tools//src/conditions:darwin": ["libarb-link.dylib"],
        "//conditions:default": ["libarb-link.so"],
    }),
    hdrs = headers,
    includes = ["."],
    deps = [
        "//third_party/flint",
    ],
)
