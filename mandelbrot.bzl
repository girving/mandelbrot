# Build extensions

load("@local_cuda//:defs.bzl", "if_local_cuda")

copts = ["-std=c++17", "-Wall", "-Werror", "-Wsign-compare"]
nvopts = ["-std=c++17", "-Wall", "-Werror=all-warnings", "-Wsign-compare"]

def cc_tests(names, deps, data=[], size="medium", copts=copts):
  deps = deps + ["@com_google_googletest//:gtest", "@com_google_googletest//:gtest_main"]
  linkopts = ["-Wno-unused-command-line-argument"]
  for name in names:
    native.cc_test(name=name, srcs=[name + ".cc"], copts=copts, linkopts=linkopts, deps=deps,
                   data=data, size=size)

def cuda_tests(names, deps, data=[], size="medium", copts=nvopts, features=[]):
  deps = deps + ["@com_google_googletest//:gtest", "@com_google_googletest//:gtest_main"]
  deps = deps + ["@rules_cuda//cuda:cuda_runtime"]
  features = features + ["cuda", "-use_header_modules"]
  linkopts = ["-L/usr/local/cuda/lib64", "-lcudart", "-Wno-unused-command-line-argument"]
  # Borrowed from https://github.com/tensorflow/runtime/blob/master/third_party/rules_cuda/cuda/defs.bzl
  special = if_local_cuda(
      select({
          "@rules_cuda//cuda:cuda_toolchain_detected": [],
          "//conditions:default": ["@rules_cuda//cuda:unsupported_cuda_toolchain_error"],
      }),
      ["@rules_cuda//cuda:no_cuda_toolkit_error"],
  )
  for name in names:
    native.cc_test(name=name, srcs=[name + ".cc"] + special, copts=copts, linkopts=linkopts, deps=deps,
                   data=data, size=size, features=features)
