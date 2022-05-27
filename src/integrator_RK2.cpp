#include "integrator.hpp"
#include "integrator_RK2.hpp"

#include <iostream>

using namespace integrators;

integrator_RK2::integrator_RK2(int verbosity) : integrator_base(verbosity) {}

double integrator_RK2::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_waves,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_particles,
		double time, double del_t) {

	// Store previous step
	utils::Grid<double> waves_old = waves;
	utils::Grid<double> particles_old = particles;

	// first Runge-Kutta step

	// get RHS
	utils::Grid<double> vec_rhs_waves = rhs_waves(waves,particles);
	utils::Grid<double> vec_rhs_particles = rhs_particles(waves,particles);
	// Do full step
	vec_rhs_waves *= del_t;
	vec_rhs_particles *= del_t;
	// Apply changes
	waves += vec_rhs_waves;
	particles += vec_rhs_particles;


	// second Runge-Kutta step
	// get RHS
	vec_rhs_waves = rhs_waves(waves,particles);
	vec_rhs_particles = rhs_particles(waves,particles);
	// Do half step
	vec_rhs_waves *= 0.5*del_t;
	vec_rhs_particles *= 0.5*del_t;

	waves_old *= 0.5;
	particles_old *= 0.5;

	waves *= 0.5;
	particles *= 0.5;

	waves += waves_old;
	waves += vec_rhs_waves;
	particles += particles_old;
	particles += vec_rhs_particles;

//	waves += del_t*vec_rhs_waves;
//	particles += del_t*vec_rhs_particles;

	time += del_t;
	return time;
}
