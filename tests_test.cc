// Tests of the test infrastructure

#include "tests.h"
namespace mandelbrot {
namespace {

TEST(assert_in_if) {
  // ASSERTs are single expressions, so they nest in unbraced if/else without dangling-else problems
  const bool yes = true, no = false;
  if (yes)
    ASSERT_EQ(1, 1) << "not printed";
  else
    ASSERT_TRUE(false);
  if (no)
    ASSERT_TRUE(false);
  // Messages are only evaluated on failure
  int evaluated = 0;
  const auto count = [&]() { evaluated++; return "message"; };
  ASSERT_EQ(1, 1) << count();
  ASSERT_TRUE(true) << count();
  ASSERT_EQ(evaluated, 0);
  ASSERT_THROW(ASSERT_EQ(1, 2) << count(), test_error);
  ASSERT_EQ(evaluated, 1);
  // Failures still throw, with or without a streamed message
  ASSERT_THROW(ASSERT_EQ(1, 2), test_error);
  ASSERT_THROW(ASSERT_LT(2, 1) << " (expected failure message)", test_error);
  ASSERT_THROW({ if (yes) ASSERT_FALSE(true); }, test_error);
}

}  // namespace
}  // namespace mandelbrot
