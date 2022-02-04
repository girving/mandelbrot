// Preprocessor utilities
#pragma once

#define UNPAREN(...) __VA_ARGS__
#define COMMA_UNPAREN(...) __VA_OPT__(,) __VA_ARGS__
