#include "galaxy.hpp"
#include "params.hpp"

int main() {
  auto params = utils::Params("");
  auto galaxy = core::Galaxy(params);
  galaxy.build_space_axis();
  galaxy.build_wavenumber_axis();
  galaxy.build_momentum_axis();
  galaxy.build_initial_condition();
  galaxy.build_wave_source();
  galaxy.build_cr_source();
  galaxy.build_advection();
  galaxy.dump_profiles();
}