// Print formatted content with a newline
#pragma once

#include "format.h"
namespace mandelbrot {

using std::ostringstream;

// Also print output to given file
void tee(const string& path);

void print();
void print(const string& s);
void print_error(const string& s);

template<class T> string str(const T& x) { ostringstream s; s << x; return s.str(); }

template<class T> static inline void print(const T& x) { print(str(x)); }
template<class... Args> void print(const char* fmt, const Args&... args) { print(format(fmt, args...)); }

template<class T> static inline void print_error(const T& x) { print_error(str(x)); }
template<class... Args> void print_error(const char* fmt, const Args&... args) { print_error(format(fmt, args...)); }

}  // namespace mandelbrot
