#include "shock_acceleration.hpp"
#include "integrator_euler.hpp"
#include "integrator_RK2.hpp"
#include "integrator_RK4.hpp"

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


	std::cout << " Building an integrator " << "\n";

	TheIntegrator my_integrator = TheIntegrator::RK2;
	int verbosity = 10;

	std::unique_ptr<integrators::integrator_base> integrator;

	if(my_integrator == TheIntegrator::Euler) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_euler(verbosity));
	} else if(my_integrator == TheIntegrator::RK2) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_RK2(verbosity));
	} else if(my_integrator == TheIntegrator::RK4) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_RK4(waves, verbosity));
	}


	std::cout << " Setting up the problem \n";

	shocks::shock_acceleration accelerator(Nx,Np);
	double time = 0.;
	double del_t = 0.01;
	int i_time = 0;

	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_waves = std::bind(&shocks::shock_acceleration::get_rhs_waves, accelerator, std::placeholders::_1, std::placeholders::_2);
	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_particles = std::bind(&shocks::shock_acceleration::get_rhs, accelerator, std::placeholders::_1, std::placeholders::_2);

//	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double &>)> rhs_waves = [=](utils::Grid<double> &a, utils::Grid<double> &b) {
//		accelerator.get_rhs( a, b);
//	}


//	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_waves = f_rhs_time_exp_waves;

	double epsilon0 = 1.e-3;
	for(int i_time=0; i_time<10; ++i_time) {
	std::cout << " particles " << particles.get(50, 2) << "\n";
	time = integrator->step(waves, particles,
					rhs_waves, rhs_particles, time, del_t);

	std::cout << " particles " << particles.get(50, 2) << "\n";
	}

	return 0;
}
