# Mandelbrot area

workspace(name = "mandelbrot")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "tinyformat",
    urls = [
        "https://github.com/c42f/tinyformat/archive/3a33bbf65442432277eee079e83d3e8fac51730c.tar.gz",
    ],
    sha256 = "52c7b9cb9558f57fbfbdbcfbb9d956793a475886a0f7a21632115978cdd7f8be",
    strip_prefix = "tinyformat-3a33bbf65442432277eee079e83d3e8fac51730c",
    build_file = "//third_party:tinyformat.BUILD",
)

http_archive(
    name = "com_google_googletest",
    urls = [
        "https://github.com/google/googletest/archive/d175c8bf823e709d570772b038757fadf63bc632.tar.gz",
    ],
    sha256 = "39a708e81cf68af02ca20cad879d1dbd055364f3ae5588a5743c919a51d7ad46",
    strip_prefix = "googletest-d175c8bf823e709d570772b038757fadf63bc632",
    build_file = "//third_party:googletest.BUILD",
)

http_archive(
    name = "arb",
    urls = [
        "https://github.com/fredrik-johansson/arb/archive/refs/tags/2.22.0.tar.gz",
    ],
    sha256 = "3e40ab8cf61c0cd63d5901064d73eaa2d04727bbdc6eebb1727997958a14f24d",
    strip_prefix = "arb-2.22.0",
    build_file = "//third_party:arb.BUILD",
)

http_archive(
    name = "flint",
    urls = [
        "https://github.com/wbhart/flint2/archive/refs/tags/v2.8.4.tar.gz",
    ],
    sha256 = "cd80a6de60bae54f68a17eac869dbd2cbf3a07b8612bbdacd83ae68ce048d60d",
    strip_prefix = "flint2-2.8.4",
    build_file = "//third_party:flint.BUILD",
)

http_archive(
    name = "mpfr",
    urls = [
        "https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.xz",
    ],
    sha256 = "0c98a3f1732ff6ca4ea690552079da9c597872d30e96ec28414ee23c95558a7f",
    strip_prefix = "mpfr-4.1.0",
    build_file = "//third_party:mpfr.BUILD",
)

#http_archive(
#    name = "gmp",
#    urls = [
#        "https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz",
#        "https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz",
#    ],
#    sha256 = "fd4829912cddd12f84181c3451cc752be224643e87fac497b69edddadc49b4f2",
#    strip_prefix = "gmp-6.2.1",
#    build_file = "//third_party:gmp.BUILD",
#)
