// Declare a value lives on the GPU
#pragma once

#include <type_traits>
#include <utility>
namespace mandelbrot {

using std::conditional_t;
using std::is_const_v;
template<class S> struct Complex;
template<class T> struct IsDeviceT;
template<class T> struct UndeviceT;
template<class T> struct AddComplexT;

// Device<T> is a T that lives on the GPU
template<class T_> struct Device {
  typedef T_ T;
private:
  T x;
};

// Type traits
template<class T> static constexpr bool is_device = IsDeviceT<T>::value;
template<class T> using Undevice = typename UndeviceT<T>::type;
template<class T> using AddComplex = typename AddComplexT<T>::type;

// Type trait implementations
template<class T> struct IsDeviceT { static constexpr bool value = false; };
template<class T> struct IsDeviceT<Device<T>> { static constexpr bool value = true; };
template<class T> struct IsDeviceT<const T> : public IsDeviceT<T> {};
template<class T> struct IsDeviceT<T*> : public IsDeviceT<T> {};
template<class T> struct UndeviceT { typedef T type; };
template<class T> struct UndeviceT<Device<T>> { typedef T type; };
template<class T> struct AddComplexT { typedef Complex<T> type; };
template<class T> struct AddComplexT<const T> { typedef const AddComplex<T> type; };
template<class T> struct AddComplexT<Device<T>> { typedef Device<Complex<T>> type; };

// Change Device<T>* to T*
template<class T> static inline auto undevice(T&& x) { return std::forward<T>(x); }  // Do nothing by default
template<class T> static inline T* undevice(Device<T>* p) { return reinterpret_cast<T*>(p); }
template<class T> static inline const T* undevice(const Device<T>* p) { return reinterpret_cast<const T*>(p); }

}  // namespace mandelbrot
