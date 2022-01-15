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
