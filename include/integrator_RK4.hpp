#ifndef INTEGRATOR_RK4_HPP
#define INTEGRATOR_RK4_HPP

#include "integrator.hpp"

namespace integrators {
class integrator_RK4 : public integrator_base {
public:
	integrator_RK4(utils::Grid<double> &waves, int verbosity);

	double step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
			std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
			std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
			double time, double del_t);
private:
	utils::Grid<double> k1_waves, k2_waves, k3_waves, k4_waves;
	utils::Grid<double> k1_particles, k2_particles, k3_particles, k4_particles;
	utils::Grid<double> waves_star, particles_star;

};
}  // namespace integrators

#endif
