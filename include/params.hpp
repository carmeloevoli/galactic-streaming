#ifndef PARAMS_H
#define PARAMS_H

#include <utility>

#include "units.hpp"

using Range = std::pair<double, double>;

namespace utils {

class Params {
 private:
  Range m_rigidity_range = {1e17 * SI::eV, 1e21 * SI::eV};
  double m_halo_size = 4. * SI::kpc;

 public:
  explicit Params(const char* inputFilename);
  virtual ~Params() = default;
  void print();

  const Range& rigidity_range = m_rigidity_range;
  const double& halo_size = m_halo_size;
};

}  // namespace utils

#endif  // SIMPROP_PARAMS_H