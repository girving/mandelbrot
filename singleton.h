// Singletons that gracefully shutdown
#pragma once

#include "noncopyable.h"
#include "shutdown.h"
#include <type_traits>
namespace mandelbrot {

using std::is_base_of_v;

// Make sure we don't accidentally create a non-singleton
class Single {
  template<class D> friend struct Singleton;
  Single() = default;
};

template<class Derived> struct Singleton : public Noncopyable {
protected:
  Singleton(Single) {
    // Call Derived::clear() on shutdown
    on_shutdown([]() { single().clear(); });
  }
public:

  static Derived& single() {
    static_assert(is_base_of_v<Singleton, Derived>);
    static Derived single((Single()));
    return single;
  }
};

}  // namespace mandelbrot
