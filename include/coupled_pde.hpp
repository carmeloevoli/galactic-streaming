#ifndef COUPLED_PDE
#define COUPLED_PDE

#include "grid.hpp"

/**
 * Virtual base class for coupled PDEs
 * Derived classes provide a right-hand side for both components and a routine for
 * boundary conditions. For the latter, the base class provides a default version
 * */
namespace coupled_pdes {

class coupled_pde {
public:
	coupled_pde() {}
	virtual ~coupled_pde() {}
	/**
	 * right-hand side for first coupled PDE component, i.e. waves
	 * */
	virtual utils::Grid<double> get_rhs_waves(utils::Grid<double> &waves,
			utils::Grid<double> &particles) = 0;

	/**
	 * right-hand side for second coupled PDE component, i.e. particles
	 * */
	virtual utils::Grid<double> get_rhs_particles(utils::Grid<double> &waves,
			utils::Grid<double> &particles) = 0;

	/**
	 * Boundary conditions for both components of coupled PDE
	 * Default version with zero-flux and zero-wavepower boundaries
	 * To be overloaded
	 * */
	virtual void set_BCs(utils::Grid<double> &waves,
			utils::Grid<double> &particles);
};

}


#endif
