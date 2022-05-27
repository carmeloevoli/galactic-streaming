#ifndef GALAXY_H
#define GALAXY_H

#include <vector>

#include "grid.hpp"
#include "params.hpp"

namespace core {

class Galaxy {
 protected:
  utils::Params m_params;

  size_t m_Nz;
  size_t m_Nk;

  std::vector<double> m_R;
  std::vector<double> m_k;
  std::vector<double> m_z;
  std::vector<double> m_vA;

  utils::Grid<double> m_W;
  utils::Grid<double> m_f;
  utils::Grid<double> m_Q_W;
  utils::Grid<double> m_Q_f;

  //   TGrid2D<double> Q_w;
  //   TGrid2D<double> Q_cr;

  //   // TGrid2D<double> W_ext;
  //   TGrid2D<double> D_kk;
  //   TGrid2D<double> D_zz;
  //   TGrid2D<double> dfdz;
  //   TGrid2D<double> v_A;
  //   TGrid2D<double> dpdt;
  //   TGrid2D<double> B_z;
  //   TGrid2D<double> Gamma_CR;
  //   TGrid2D<double> WGamma_CR;
  //   TGrid2D<double> fcr;

 public:
  Galaxy(const utils::Params& params) : m_params(params) {}

  void build_momentum_axis();
  void build_rigidity_axis();
  void build_space_axis();
  void build_initial_condition();
  void build_wave_source();
  void build_cr_source();
  void build_advection();
  //   void build_energy_losses();
  //   void build_magnetic_field();
};

}  // namespace core

#endif