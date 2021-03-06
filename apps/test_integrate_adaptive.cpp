#include "integrators/integrator_euler.hpp"
#include "integrators/integrator_RK2.hpp"
#include "integrators/integrator_RK4.hpp"

#include <iostream>
#include <memory>
#include <math.h>

using namespace integrators;

class time_adaptive_test : public coupled_pdes::coupled_pde {
public:
	time_adaptive_test() {}

utils::Grid<double> get_rhs_waves(utils::Grid<double> &waves,utils::Grid<double> &particles) {
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	double lambd = -2.;
	utils::Grid<double> rhs_w(Nx, Nz);
	rhs_w = waves*lambd;
	return rhs_w;
}

utils::Grid<double> get_rhs_particles(utils::Grid<double> &waves,utils::Grid<double> &particles) {
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	double lambd = -2.;
	utils::Grid<double> rhs_p(Nx, Nz);
	rhs_p = particles*lambd;
	return rhs_p;
}

};


enum class TheIntegrator {
	Euler,
	RK2,
	RK4
};

int main() {

	std::shared_ptr<time_adaptive_test> adaptive = std::make_shared<time_adaptive_test>();


	utils::Grid<double> waves(3,3,1.); // set value to 1
	utils::Grid<double> particles(3,3,2.); // set values to 2

	TheIntegrator my_integrator = TheIntegrator::RK2;
	int verbosity = 10;

	std::unique_ptr<integrators::integrator_base> integrator;

	if(my_integrator == TheIntegrator::Euler) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_euler(adaptive, verbosity));
	} else if(my_integrator == TheIntegrator::RK2) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_RK2(adaptive, verbosity));
	} else if(my_integrator == TheIntegrator::RK4) {
		integrator = std::unique_ptr< integrators::integrator_base > (new integrator_RK4(adaptive, waves, verbosity));
	}

	double time = 0.;
	double del_t = 0.01;
	int i_time = 0;

	double epsilon0 = 1.e-9;
	while(time<1.) {

		time = integrator->step_adaptive(waves, particles, time, del_t, epsilon0);
		if(i_time%1==0) {
			std::cout << " step " << i_time << " at " << time << " ";
			std::cout << " with dt: " << del_t << " -> ";
			std::cout << waves.get_value(1, 1) << " ";
			std::cout << waves.get_value(1, 1) - exp(-2.*time);
			std::cout << "\n";
			std::cout << "                               -> ";
			std::cout << particles.get_value(1, 1) << " ";
			std::cout << 2.*exp(-2.*time) << " ";
			std::cout << particles.get_value(1, 1) - 2.*exp(-2.*time);
			std::cout << "\n";
		}
//		if(i_time==2) exit(2);
		i_time++;

	}
	std::cout << " Finished after " << i_time << " steps\n";
	return 0;
}
