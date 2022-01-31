// For assertions that scalars aren't intervals
#pragma once

namespace mandelbrot {

template<int n> struct Expansion;

template<class T> struct IsIntervalT;
template<class T> struct IsIntervalT<const T> : public IsIntervalT<T> {};
template<> struct IsIntervalT<double> { static constexpr bool value = false; };
template<int n> struct IsIntervalT<Expansion<n>> { static constexpr bool value = false; };

template<class T> static constexpr bool is_interval = IsIntervalT<T>::value;

}  // namespace mandelbrot
