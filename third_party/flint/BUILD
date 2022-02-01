# GMP

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # LGPL

headers = [
    "flint/config.h",
    "flint/flint.h",
    "flint/flint-config.h",
    "flint/d_mat.h",
    "flint/d_vec.h",
    "flint/double_extras.h",
    "flint/exception.h",
    "flint/fmpq.h",
    "flint/fmpq_poly.h",
    "flint/fmpz.h",
    "flint/fmpz-conversions.h",
    "flint/fmpz_factor.h",
    "flint/fmpz_mat.h",
    "flint/fmpz_poly.h",
    "flint/fmpz_poly_factor.h",
    "flint/fmpz_vec.h",
    "flint/gmpcompat.h",
    "flint/longlong.h",
    "flint/mpf_mat.h",
    "flint/mpf_vec.h",
    "flint/nmod_mat.h",
    "flint/nmod_poly.h",
    "flint/nmod_poly_factor.h",
    "flint/nmod_vec.h",
    "flint/thread_pool.h",
    "flint/thread_support.h",
    "flint/ulong_extras.h",
]

genrule(
    name = "headers",
    outs = headers,
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "\n".join([
            "for h in $(OUTS); do",
            "  s=`brew --prefix`/include/flint/$$(basename $$h)",
            "  if [ -e $$s ]; then",
            "    ln -s `brew --prefix`/include/flint/$$(basename $$h) $$h",
            "  else",
            "    echo nonexistent > $$s",
            "  fi",
            "done",
        ]),
        "//conditions:default": "\n".join([
            "if [ -e /usr/include/flint ]; then",
            "  src=/usr/include",
            "elif [ -e /usr/include/x86_64-linux-gnu/flint ]; then",
            "  src=/usr/include/x86_64-linux-gnu",
            "fi",
            "for h in $(OUTS); do",
            "  s=$$src/flint/$$(basename $$h)",
            "  if [ -e $$s ]; then",
            "    ln -s $$src/flint/$$(basename $$h) $$h",
            "  else",
            "    touch $$h",
            "  fi",
            "done",
        ]),
    }),
)

genrule(
    name = "dylib",
    outs = ["libflint-link.dylib"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "ln -s `brew --prefix`/lib/libflint.dylib $(OUTS)",
        "//conditions:default": "touch $(OUTS)",
    }),
)

genrule(
    name = "so",
    outs = ["libflint-link.so"],
    cmd = select({
        "@bazel_tools//src/conditions:darwin": "touch $(OUTS)",
        "//conditions:default": "\n".join([
            "if [ -e /usr/lib/libflint.so ]; then",
            "  ln -s /usr/lib/libflint.so $(OUTS)",
            "elif [ -e /usr/lib/x86_64-linux-gnu/libflint.so ]; then",
            "  ln -s /usr/lib/x86_64-linux-gnu/libflint.so $(OUTS)",
            "fi",
        ]),
    }),
)

cc_library(
    name = "flint",
    srcs = select({
        "@bazel_tools//src/conditions:darwin": ["libflint-link.dylib"],
        "//conditions:default": ["libflint-link.so"],
    }),
    hdrs = headers,
    includes = ["."],
    deps = [
        "//third_party/gmp",
        "//third_party/mpfr",
    ],
)
