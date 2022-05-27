#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_integration.h>

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

inline double momentum_2_radius(double R, double B) {
  return R / SI::elementaryCharge / SI::cLight / B;
}

inline double radius_2_momentum(double r_L, double B) {
  return SI::cLight * SI::elementaryCharge * B * r_L;
}

inline double magnetic_energy_density(double B) {
  const auto mu_0 = 1.25663706212e-6;
  return .5 * B * B / mu_0;
}

inline double Gaussian(double x, double sigma) {
  const auto sigma2 = sigma * sigma;
  return std::pow(2. * M_PI * sigma2, -0.5) * std::exp(-(x * x) / 2. / sigma2);
}

inline double power_law_with_cutoff(double x, double slope, double x_cutoff) {
  return std::pow(x, -slope) * std::exp(-x / x_cutoff);
}

double f_source_integral(double x, void* params) {
  double alpha = *(double*)params;
  double f = pow(x, 2. - alpha) * (sqrt(x * x + 1.) - 1.);
  return f;
}

double source_integral(double alpha, double x_max) {
  int LIMIT = 10000;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  double result, error;
  gsl_function F;
  F.function = &f_source_integral;
  F.params = &alpha;
  gsl_integration_qag(&F, 1, x_max, 0, 1e-6, LIMIT, 3, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

template <typename T>
size_t closest(std::vector<T> const& vec, T value) {
  auto const it = std::lower_bound(vec.begin(), vec.end(), value);
  if (it == vec.end()) {
    return -1;
  }
  return it - vec.begin();
}

}  // namespace utils

#endif