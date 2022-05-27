#ifndef PARAMS_H
#define PARAMS_H

#include "units.hpp"

namespace utils {

class Params {
 private:
  double m_halo_size = 4. * SI::kpc;
  double m_k_min = 1e-2 / SI::pc;
  double m_k_max = 1e7 / SI::pc;
  size_t m_z_size = 80;
  size_t m_k_size = 9 * 32;
  double m_magnetic_field = 3. * SI::muG;
  double m_eta_B = 0.1;
  double m_k_0 = 1e-1 / SI::pc;

 public:
  explicit Params(const char* inputFilename);
  virtual ~Params() = default;
  void print();

  const double& halo_size = m_halo_size;
  const double& k_min = m_k_min;
  const double& k_max = m_k_max;
  const size_t& z_size = m_z_size;
  const size_t& k_size = m_k_size;
  const double& magnetic_field = m_magnetic_field;
  const double& k_0 = m_k_0;
  const double& eta_B = m_eta_B;
};

}  // namespace utils

#endif  // SIMPROP_PARAMS_H