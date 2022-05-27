#include "galaxy.hpp"

#include <cassert>
#include <iostream>

#include "utils.hpp"

namespace core {

void Galaxy::build_space_axis() {
  m_z = utils::build_lin_axis<double>(-m_params.halo_height, m_params.halo_height, m_params.z_size);
  m_Nz = m_params.z_size;
  std::cout << "built z axis between " << m_z.front() / SI::kpc << " kpc and "
            << m_z.back() / SI::kpc << " kpc with " << m_Nz << " points.\n";
}

void Galaxy::build_momentum_axis() {
  m_k = utils::build_log_axis<double>(m_params.k_min, m_params.k_max, m_params.k_size);
  m_Nk = m_params.k_size;
  std::cout << "built k axis between " << m_k.front() * SI::pc << " pc^-1 and "
            << m_k.back() * SI::pc << " pc^-1 with " << m_Nk << " points.\n";
}

void Galaxy::build_rigidity_axis() {
  assert(m_k.size() == m_Nk);
  auto min_rigidity = utils::radius_2_rigidity(1. / m_params.k_max, m_params.magnetic_field);
  auto max_rigidity = utils::radius_2_rigidity(1. / m_params.k_min, m_params.magnetic_field);
  m_R = utils::build_log_axis<double>(min_rigidity, max_rigidity, m_Nk);
  std::cout << "built R axis between " << m_R.front() / SI::GV << " GV and " << m_R.back() / SI::GV
            << " GV with " << m_Nk << " points.\n";
}

void Galaxy::build_initial_condition() {
  const double TINY = 1e-4;
  // build initial W with a Kolmorogov spectrum
  m_W = utils::Grid<double>(m_Nk, m_Nz);
  for (size_t ik = 0; ik < m_Nk; ++ik) {
    auto value = 2.0 * m_params.eta_B / 3.0 / m_params.k_0;
    if (m_k.at(ik) > m_params.k_0) value *= std::pow(m_k.at(ik) / m_params.k_0, -5. / 3.);
    for (size_t iz = 0; iz < m_Nz; ++iz) {
      m_W.get(ik, iz) = TINY * value;
    }
  }
  // build initial empty f
  m_f = utils::Grid<double>(m_Nk, m_Nz, 0.);
}

void Galaxy::build_wave_source() {
  m_Q_W = utils::Grid<double>(m_Nk, m_Nz);
  auto U_B = utils::magnetic_energy_density(m_params.magnetic_field);
  auto disk_area = M_PI * pow2(m_params.disk_radius);
  auto luminosity = m_params.eta_wave * m_params.E_SN * m_params.R_SN;
  auto q_0 = luminosity / U_B / disk_area;
  auto h = m_params.disk_height;
  auto sigma_k = 0.1 * m_params.k_0;
  for (size_t i = 0; i < m_Nk; ++i) {
    auto gauss_k = utils::Gaussian(m_k[i] - m_params.k_0, sigma_k);
    for (size_t j = 0; j < m_Nz; ++j) {
      auto gauss_z = utils::Gaussian(m_z[j], h);
      m_Q_W.get(i, j) = q_0 * gauss_k * gauss_z;
    }
  }
}

// double compute_constant_CR_source_term() {
// 	double out = params.efficiency_cr_sn() * params.kinetic_energy_sn() * params.rate_sn();
// 	out /= 4. * M_PI * params.source_disk_surface() * I_of_alpha(params.alpha()) * c_light *
// pow4(mass_proton_c); 	return out;
// }

void Galaxy::build_cr_source() {
  m_Q_f = utils::Grid<double>(m_Nk, m_Nz);
  const double q_0 = 1;  // compute_constant_CR_source_term();
  const auto h = m_params.disk_height;
  for (size_t i = 0; i < m_Nk; ++i) {
    auto x = m_R[i] / m_params.R_0;
    auto spectrum = utils::power_law_with_cutoff(x, m_params.slope, m_params.source_cutoff);
    for (size_t j = 0; j < m_Nz; ++j) {
      auto gauss_z = utils::Gaussian(m_z[j], h);
      m_Q_f.get(i, j) = q_0 * spectrum * gauss_z;
    }
  }
}

void Galaxy::build_advection() {
  const auto h = 25.0 * SI::pc;
  m_vA.reserve(m_Nz);
  for (size_t j = 0; j < m_Nz; ++j) {
    m_vA.emplace_back(2.0 * atan(m_z[j] / h) / M_PI * m_params.vA_infty);
  }
}

}  // namespace core