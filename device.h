// Declare a value lives on the GPU
#pragma once

namespace mandelbrot {

// Device<T> is a T that lives on the GPU
template<class T> struct Device {
private:
  T x;
};

template<class T> struct IsDeviceT { static constexpr bool value = false; };
template<class T> struct IsDeviceT<Device<T>> { static constexpr bool value = true; };
template<class T> struct IsDeviceT<const Device<T>> { static constexpr bool value = true; };
template<class T> static constexpr bool is_device = IsDeviceT<T>::value;

template<class T> struct UndeviceT { typedef T type; };
template<class T> struct UndeviceT<Device<T>> { typedef T type; };
template<class T> using Undevice = typename UndeviceT<T>::type;

}  // namespace mandelbrot
