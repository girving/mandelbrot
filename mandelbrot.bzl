# Build extensions

copts = ["-std=c++2a", "-Wall", "-Werror", "-Wsign-compare"]

def cc_tests(names, deps, data=[], size="medium", copts=copts):
  deps = deps + ["@com_google_googletest//:gtest", "@com_google_googletest//:gtest_main"]
  linkopts = ["-Wno-unused-command-line-argument"]
  for name in names:
    native.cc_test(name=name, srcs=[name + ".cc"], copts=copts, linkopts=linkopts, deps=deps,
                   data=data, size=size)
