#include "coupled_pde.hpp"

namespace coupled_pdes {

void coupled_pde::set_BCs(utils::Grid<double> &waves,
			utils::Grid<double> &particles) {
	// Default version assumes zero flux and zero wave power at all boundaries

	// Boundaries for waves
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	// Lower and upper x-boundary
	for(size_t i_z=0; i_z<Nz; ++i_z) {
		waves.set_value(0   , i_z, 0.);
		waves.set_value(Nx-1, i_z, 0.);
	}
	// Lower and upper z-boundary
	for(size_t i_x=0; i_x<Nx; ++i_x) {
		waves.set_value(i_x,0   , 0.);
		waves.set_value(i_x,Nz-1, 0.);
	}


	// Boundaries for particles
	Nx=particles.get_Nx();
	Nz=particles.get_Nz();
	// Lower and upper x-boundary
	for(size_t i_z=0; i_z<Nz; ++i_z) {
		particles.set_value(0   , i_z, 0.);
		particles.set_value(Nx-1, i_z, 0.);
	}
	// Lower and upper z-boundary
	for(size_t i_x=0; i_x<Nx; ++i_x) {
		particles.set_value(i_x,0   , 0.);
		particles.set_value(i_x,Nz-1, 0.);
	}
}

}
