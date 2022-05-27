#ifndef UTILS_H
#define UTILS_H

#include <stdexcept>

#include "units.hpp"

namespace utils {

template <typename T>
std::vector<T> build_lin_axis(const T& min, const T& max, const size_t& size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const auto dx = (max - min) / (T)(size - 1);
  std::vector<T> v(size);
  for (size_t i = 0; i < size; ++i) {
    const auto value = min + dx * i;
    v[i] = value;
  }
  return v;
}

template <typename T>
std::vector<T> build_log_axis(const T& min, const T& max, const size_t& size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const auto delta_log = std::exp(std::log(max / min) / (T)(size - 1));
  std::vector<T> v(size);
  for (size_t i = 0; i < size; ++i) {
    const auto value = std::exp(std::log(min) + (T)i * std::log(delta_log));
    v[i] = value;
  }
  return v;
}

inline double larmor_radius(double R, double B) { return R / SI::cLight / B; }

inline double inverse_larmor_radius(double k, double B) { return SI::cLight * B / k; }

}  // namespace utils

#endif