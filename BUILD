package(default_visibility = ["//visibility:public"])

cc_binary(
    name = "area",
    srcs = [
      "arb.h",
      "arb.cc",
      "arf.h",
      "arf.cc",
      "area.cc",
      "poly.h",
      "poly.cc",
      "print.h",
      "wall_time.h",
    ],
    copts = ["-std=c++20", "-Wall", "-Werror"],
    deps = [
        ":known",
        "//third_party/arb",
        "@tinyformat//:tinyformat",
    ],
)

# DO NOT SUBMIT
#cc_binary(
#    name = "area-c",
#    srcs = ["area.c"],
#    copts = ["-std=c99", "-Wall", "-Werror"],
#    deps = [
#        ":known",
#        "//third_party/arb",
#    ],
#)

cc_library(
    name = "known",
    hdrs = ["known.h"],
    copts = ["-Wall", "-Werror"],
)
