// Mandelbrot area via custom power series

#include "series.h"
#include <tuple>
namespace mandelbrot {

using std::tuple;

// Estimate area given computed f (warning: not g)
template<class S> S area(SeriesView<S> f);

// f = 1, so g = log f = 0 + O(z)
template<class T> Series<T> bottcher_base();

// Double the number of terms of g, assuming existing terms are correct.
// Returns f = exp g and the resulting error estimate.
template<class T> tuple<Series<T>,Undevice<T>> bottcher_step(Series<T>& g, const double tol);

// Compute series up to 1<<max_k terms, printing and checking errors along the way
template<class T> void areas(const int max_k, const double tol);

}  // namespace mandelbrot
