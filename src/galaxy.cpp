#include "galaxy.hpp"

#include <cassert>
#include <iostream>

#include "utils.hpp"

namespace core {

void Galaxy::build_space_axis() {
  m_z = utils::build_lin_axis<double>(-m_params.halo_size, m_params.halo_size, m_params.z_size);
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
  const auto min_rigidity = utils::inverse_larmor_radius(m_params.k_max, m_params.magnetic_field);
  const auto max_rigidity = utils::inverse_larmor_radius(m_params.k_min, m_params.magnetic_field);
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

}  // namespace core