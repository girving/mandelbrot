// Very simple testing infrastructure
#pragma once

#include "debug.h"
#include "format.h"
#include "noncopyable.h"
#include "print.h"
namespace mandelbrot {

using std::cout;
using std::exception;
using std::function;
using std::ostream;
using std::string_view;

int register_test(const char* name, void (*test)());
void test_throw(const function<void()>& f, const char* sx, const char* se, const char* function, const int line);
struct test_error : public exception {};
string red(const string& s);
string green(const string& s);
string blue(const string& s);

#define TEST(name) \
  static void name##_test(); \
  __attribute__((unused)) const int name##_ignored = register_test(#name, name##_test); \
  void name##_test()

// On failure, ASSERTs stream an optional message into cout, then throw.  Each ASSERT is a single expression
//   ({ check; }) ? (void)0 : DieOnAssign() = cout << message...
// so it nests safely in unbraced if/else, and since the trailing << binds inside the third operand of ?:, the
// message is only evaluated on failure.
struct DieOnAssign {
  void operator=(const ostream&) { print(); throw test_error(); }
};

template<class X> static inline bool
test_bool(const char* sx, X&& x, const bool y, const char* function, const int line) {
  const bool r = bool(x) == y;
  if (!r)
    print("%s %s != %s", red(tfm::format("%s:%d:", function, line)), sx, y ? "true" : "false");
  return r;
}

template<class X,class Y,class Op> static inline bool
test_compare(const char* sx, X&& x, const char* sy, Y&& y, Op&& op, const char* nop,
             const char* function, const int line) {
  const bool r = op(x, y);
  if (!r)
    print("%s %s %s %s (%s %s %s)", red(tfm::format("%s:%d:", function, line)), sx, nop, sy, x, nop, y);
  return r;
}

#define ASSERT2(x, y, op, nop) \
  ({ test_compare(#x, (x), #y, (y), op, nop, __FUNCTION__, __LINE__); }) ? (void)0 : DieOnAssign() = cout

#define ASSERT_TRUE(x) \
  ({ test_bool(#x, (x), true, __FUNCTION__, __LINE__); }) ? (void)0 : DieOnAssign() = cout
#define ASSERT_FALSE(x) \
  ({ test_bool(#x, (x), false, __FUNCTION__, __LINE__); }) ? (void)0 : DieOnAssign() = cout

#define ASSERT_EQ(x, y) ASSERT2(x, y, std::equal_to(), "!=")
#define ASSERT_NE(x, y) ASSERT2(x, y, std::not_equal_to(), "!!=")
#define ASSERT_LE(x, y) ASSERT2(x, y, std::less_equal(), "!<=")
#define ASSERT_LT(x, y) ASSERT2(x, y, std::less(), "!<")

#define ASSERT_THROW(x, e) \
  test_throw([&]() { x; }, #x, #e, __FUNCTION__, __LINE__)

struct Tmpfile : public Noncopyable {
  const string path;
  Tmpfile(string_view prefix);
  ~Tmpfile();
};

struct Tmpdir : public Noncopyable {
  const string path;
  Tmpdir(string_view prefix);
  ~Tmpdir();
};

}  // namespace mandelbrot
