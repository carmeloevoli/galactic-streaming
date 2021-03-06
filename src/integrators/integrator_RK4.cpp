#include "integrators/integrator.hpp"
#include "integrators/integrator_RK4.hpp"

#include <iostream>

using namespace integrators;

integrator_RK4::integrator_RK4(std::shared_ptr<coupled_pdes::coupled_pde> pde,
		utils::Grid<double> &waves, int verbosity):
		integrator_base(pde,verbosity) {
	scheme_order = 4.;

	int Nx = waves.get_Nx();
	int Nz = waves.get_Nz();

	k1_waves.set_grid_size(Nx, Nz);
	k2_waves.set_grid_size(Nx, Nz);
	k3_waves.set_grid_size(Nx, Nz);
	k4_waves.set_grid_size(Nx, Nz);

	k1_particles.set_grid_size(Nx, Nz);
	k2_particles.set_grid_size(Nx, Nz);
	k3_particles.set_grid_size(Nx, Nz);
	k4_particles.set_grid_size(Nx, Nz);

	waves_star.set_grid_size(Nx, Nz);
	particles_star.set_grid_size(Nx, Nz);
}

double integrator_RK4::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		double time, double del_t) {

	double sixth = 1./6.;

	// Store previous step
	utils::Grid<double> waves_old = waves;
	utils::Grid<double> particles_old = particles;

	// first Runge-Kutta step

	// get RHS
	utils::Grid<double> vec_rhs_waves = m_pde->get_rhs_waves(waves,particles);
	utils::Grid<double> vec_rhs_particles = m_pde->get_rhs_particles(waves,particles);

	// Do first step
	k1_waves = m_pde->get_rhs_waves(waves,particles);
	k1_particles = m_pde->get_rhs_particles(waves,particles);

	waves_star= waves + k1_waves*0.5*del_t;
	particles_star= particles + k1_particles*0.5*del_t;
	// Apply boundary conditions:
	m_pde->set_BCs(waves_star, particles_star);

	// Do second step
	k2_waves = m_pde->get_rhs_waves(waves_star,particles_star);
	k2_particles = m_pde->get_rhs_particles(waves_star,particles_star);

	waves_star= waves + k2_waves*0.5*del_t;
	particles_star= particles + k2_particles*0.5*del_t;
	// Apply boundary conditions:
	m_pde->set_BCs(waves_star, particles_star);

	// Do third step
	k3_waves = m_pde->get_rhs_waves(waves_star,particles_star);
	k3_particles = m_pde->get_rhs_particles(waves_star,particles_star);

	waves_star= waves + k3_waves*del_t;
	particles_star= particles + k3_particles*del_t;
	// Apply boundary conditions:
	m_pde->set_BCs(waves_star, particles_star);

	// Do fourth step
	k4_waves = m_pde->get_rhs_waves(waves_star,particles_star);
	k4_particles = m_pde->get_rhs_particles(waves_star,particles_star);

	// Get solution after full step
	waves += (k1_waves + k2_waves*2. + k3_waves*2 + k4_waves)*del_t*sixth;
	particles += (k1_particles + k2_particles*2 + k3_particles*2 + k4_particles)*del_t*sixth;
	// Apply boundary conditions:
	m_pde->set_BCs(waves, particles);

	time += del_t;
	return time;
}
