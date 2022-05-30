#include "integrator.hpp"
#include "integrator_euler.hpp"

#include <iostream>

using namespace integrators;

integrator_euler::integrator_euler(std::shared_ptr<coupled_pdes::coupled_pde> pde,
		int verbosity) : integrator_base(pde, verbosity) {
	scheme_order = 1.;
}
//integrator_euler::integrator_euler(int verbosity) : integrator_base(verbosity) {
//	scheme_order = 1.;
//}

double integrator_euler::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		double time, double del_t) {



	utils::Grid<double> vec_rhs_waves = m_pde->get_rhs_waves(waves,particles);
	utils::Grid<double> vec_rhs_particles = m_pde->get_rhs_particles(waves,particles);

	waves += vec_rhs_waves*del_t;
	particles += vec_rhs_particles*del_t;

	// Apply boundary conditions:
	m_pde->set_BCs(waves, particles);

	time += del_t;
	return time;
}

//double integrator_euler::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
//		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_waves,
//		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_particles,
//		double time, double del_t) {
//
//
//
//	utils::Grid<double> vec_rhs_waves = rhs_waves(waves,particles);
//	utils::Grid<double> vec_rhs_particles = rhs_particles(waves,particles);
//
//	waves += vec_rhs_waves*del_t;
//	particles += vec_rhs_particles*del_t;
//
//	time += del_t;
//	return time;
//}
