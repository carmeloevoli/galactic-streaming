#ifndef PARAMS_H
#define PARAMS_H

#include "units.hpp"

namespace utils {

class Params {
 private:
  double m_halo_height = 4. * SI::kpc;
  double m_disk_height = 0.1 * SI::kpc;
  double m_disk_radius = 10. * SI::kpc;
  double m_k_min = 1e-2 / SI::pc;
  double m_k_max = 1e7 / SI::pc;
  size_t m_z_size = 80;
  size_t m_k_size = 9 * 32;
  double m_magnetic_field = 3. * SI::muG;
  double m_eta_B = 0.1;
  double m_k_0 = 1e-1 / SI::pc;
  double m_slope = 4.3;
  double m_source_cutoff = SI::PV;
  double m_vA_infty = 10. * SI::km / SI::sec;
  double m_eta_wave = 0.001;
  double m_eta_cr = 0.1;
  double m_E_SN = 1e51 * SI::erg;
  double m_R_SN = 1. / (30. * SI::year);

 public:
  explicit Params(const char* inputFilename);
  virtual ~Params() = default;
  void print();

  const double& halo_height = m_halo_height;
  const double& disk_height = m_disk_height;
  const double& disk_radius = m_disk_radius;
  const double& k_min = m_k_min;
  const double& k_max = m_k_max;
  const size_t& z_size = m_z_size;
  const size_t& k_size = m_k_size;
  const double& magnetic_field = m_magnetic_field;
  const double& k_0 = m_k_0;
  const double& eta_B = m_eta_B;
  const double& slope = m_slope;
  const double& source_cutoff = m_source_cutoff;
  const double& vA_infty = m_vA_infty;
  const double& eta_wave = m_eta_wave;
  const double& eta_cr = m_eta_cr;
  const double& E_SN = m_E_SN;
  const double& R_SN = m_R_SN;
};

}  // namespace utils

#endif  // SIMPROP_PARAMS_H