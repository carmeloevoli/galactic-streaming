#include "integrator.hpp"
#include "integrator_euler.hpp"

#include <iostream>

using namespace integrators;

integrator_euler::integrator_euler(int verbosity) : integrator_base(verbosity) {}

double integrator_euler::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_waves,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_particles,
		double time, double del_t) {



	utils::Grid<double> vec_rhs_waves = rhs_waves(waves,particles);
	utils::Grid<double> vec_rhs_particles = rhs_particles(waves,particles);

	vec_rhs_waves *= del_t;
	vec_rhs_particles *= del_t;

//	waves += del_t*vec_rhs_waves;
//	particles += del_t*vec_rhs_particles;

	time += del_t;
	return time;
}
