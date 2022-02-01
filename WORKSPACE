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

# From https://github.com/tensorflow/runtime/tree/master/third_party/rules_cuda
http_archive(
    name = "rules_cuda",
    sha256 = "f80438bee9906e9ecb1a8a4ae2365374ac1e8a283897281a2db2fb7fcf746333",
    strip_prefix = "runtime-b1c7cce21ba4661c17ac72421c6a0e2015e7bef3/third_party/rules_cuda",
    urls = ["https://github.com/tensorflow/runtime/archive/b1c7cce21ba4661c17ac72421c6a0e2015e7bef3.tar.gz"],
)
load("@rules_cuda//cuda:dependencies.bzl", "rules_cuda_dependencies")
rules_cuda_dependencies()
load("@rules_cc//cc:repositories.bzl", "rules_cc_toolchains")
rules_cc_toolchains()
