#include <integrator_euler.hpp>
#include <integrator_RK2.hpp>
#include <integrator_RK4.hpp>
#include <iostream>
#include <memory>
#include <math.h>

using namespace integrators;

utils::Grid<double> f_rhs_time_exp_waves(utils::Grid<double> &waves,utils::Grid<double> &particles) {
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	double lambd = -2.;
	utils::Grid<double> rhs_w(Nx, Nz);
	rhs_w = waves*lambd;
	return rhs_w;
}

utils::Grid<double> f_rhs_time_exp_particles(utils::Grid<double> &waves,utils::Grid<double> &particles) {
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	double lambd = -2.;
	utils::Grid<double> rhs_p(Nx, Nz);
	rhs_p = particles*lambd;
	return rhs_p;
}


enum class TheIntegrator {
	Euler,
	RK2,
	RK4
};

int main() {

	utils::Grid<double> waves(3,3,1.); // set value to 1
	utils::Grid<double> particles(3,3,2.); // set values to 2

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
//	integrators::integrator_RK2 integrator(10);
////	integrators::integrator_RK2 integrator(10);

	double time = 0.;
	double del_t = 0.01;
	int i_time = 0;

	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_waves = f_rhs_time_exp_waves;
	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_particles = f_rhs_time_exp_particles;

	double epsilon0 = 1.e-9;
	while(time<1.) {


		time = integrator->step_adaptive(waves, particles,
				rhs_waves, rhs_particles, time, del_t, epsilon0);
		if(i_time%1==0) {
			std::cout << " step " << i_time << " at " << time << " ";
			std::cout << " with dt: " << del_t << " -> ";
			std::cout << waves.get_value(2, 2) << " ";
			std::cout << waves.get_value(2, 2) - exp(-2.*time);
			std::cout << "\n";
			std::cout << "                               -> ";
			std::cout << particles.get_value(2, 2) << " ";
			std::cout << 2.*exp(-2.*time) << " ";
			std::cout << particles.get_value(2, 2) - 2.*exp(-2.*time);
			std::cout << "\n";
		}
//		if(i_time==2) exit(2);
		i_time++;

	}
	std::cout << " Finished after " << i_time << " steps\n";
	return 0;
}
