#ifndef GALAXY_H
#define GALAXY_H

class Galaxy {
 protected:
  std::vector<double> R;
  std::vector<double> k;
  std::vector<double> z;
  //   TGrid2D<double> Q_w;
  //   TGrid2D<double> Q_cr;
  //   TGrid2D<double> W;
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
  Galaxy() {}

  void build_momentum_axis(const double& k_min, const double& k_max, const size_t& k_size);
  void build_rigidity_axis(const double& pc_min, const double& pc_max, const size_t& pc_size);
  void build_space_axis(const double& halo_size, const size_t& k_size);
  //   void build_wave_source_term();
  //   void build_CR_source_term();
  //   double compute_constant_CR_source_term();
  //   double compute_constant_wave_source_term();
  //   void build_W();
  //   void build_initial_condition();
  //   void build_vA();
  //   void build_energy_losses();
  //   void build_magnetic_field();
};

#endif