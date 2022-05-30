#include "shock_acceleration.hpp"
#include "integrators/integrator_euler.hpp"
#include "integrators/integrator_RK2.hpp"
#include "integrators/integrator_RK4.hpp"
#include "hdf5_output.hpp"

#include <iostream>
#include <memory>

using namespace integrators;
using namespace std::placeholders;


enum class TheIntegrator {
	Euler,
	RK2,
	RK4
};

int main() {

	// Setting the grid

	int Nx = 101, Np = 100;
	utils::Grid<double> waves(Nx,Np);
	utils::Grid<double> particles(Nx,Np);

	std::shared_ptr<shocks::shock_acceleration> accelerator = std::make_shared<shocks::shock_acceleration>(Nx,Np);

	std::cout << " Building an integrator " << "\n";

	TheIntegrator my_integrator = TheIntegrator::RK2;
	int verbosity = 10;

	std::unique_ptr<integrators::integrator_base> integrator;

	if(my_integrator == TheIntegrator::Euler) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_euler(accelerator, verbosity));
	} else if(my_integrator == TheIntegrator::RK2) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_RK2(accelerator, verbosity));
	} else if(my_integrator == TheIntegrator::RK4) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_RK4(accelerator, waves, verbosity));
	}


	std::cout << " Setting up the problem \n";

	double time = 0.;
	double del_t = 2.e-5;


	// Prepare storage operator
	hdf5_IO::hdf5_writer storage("shock_test.h5");
	storage.write_grids(accelerator->get_spatial_grid(), accelerator->get_momentum_grid());


//	double epsilon0 = 1.e-3;

	size_t i_step = 0;
	while(i_step <1000000) {

		time = integrator->step(waves, particles, time, del_t);

		if(i_step%5000==0) {
			std::cout << " step: " << i_step << " particles " << particles.get(50, 2) << "\n";
			storage.write_timestep(waves, particles, time, i_step);
		}

		i_step++;
	}

	return 0;
}
