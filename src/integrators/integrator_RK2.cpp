#include "integrators/integrator.hpp"
#include "integrators/integrator_RK2.hpp"

#include <iostream>

using namespace integrators;

integrator_RK2::integrator_RK2(std::shared_ptr<coupled_pdes::coupled_pde> pde, int verbosity) : integrator_base(pde, verbosity) {
	scheme_order = 2.;
}

double integrator_RK2::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		double time, double del_t) {

	// Store previous step
	utils::Grid<double> waves_old = waves;
	utils::Grid<double> particles_old = particles;

	// first Runge-Kutta step

	// get RHS
	utils::Grid<double> vec_rhs_waves = m_pde->get_rhs_waves(waves,particles);
	utils::Grid<double> vec_rhs_particles = m_pde->get_rhs_particles(waves,particles);

	// Do full step
	waves += vec_rhs_waves*del_t;
	particles += vec_rhs_particles*del_t;
	// Apply boundary conditions:
	m_pde->set_BCs(waves, particles);

	// second Runge-Kutta step
	// get RHS
	vec_rhs_waves = m_pde->get_rhs_waves(waves,particles);
	vec_rhs_particles = m_pde->get_rhs_particles(waves,particles);

	// Do half step
	waves = waves_old*0.5 + waves*0.5 + vec_rhs_waves*0.5*del_t;
	particles = particles_old*0.5 + particles*0.5 + vec_rhs_particles*0.5*del_t;
	// Apply boundary conditions:
	m_pde->set_BCs(waves, particles);

	time += del_t;
	return time;
}
