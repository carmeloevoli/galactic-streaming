#include "integrator.hpp"
#include "integrator_euler.hpp"

#include <iostream>

using namespace integrators;

integrator_euler::integrator_euler(int verbosity) : integrator_base(verbosity) {}

double integrator_euler::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_waves,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_particles,
		double time, double del_t) {



	utils::Grid<double> vec_rhs_waves = rhs_waves(waves,particles);
	utils::Grid<double> vec_rhs_particles = rhs_particles(waves,particles);

	std::cout << " size " << vec_rhs_waves.get_size() << "\n";

	utils::Grid<double> thingy(waves);
	thingy = vec_rhs_waves*del_t;

	std::cout << " waves " << waves.get_value(20, 50) << "\n";
	waves += vec_rhs_waves*del_t;
	std::cout << " In Euler " << vec_rhs_waves.get_value(20, 50) << " " << del_t << " ";
	std::cout << waves.get_value(20, 50) << " ";
	std::cout << thingy.get_value(20, 50) << " ";
	std::cout << "\n";
	particles += vec_rhs_particles*del_t;

	time += del_t;
	return time;
}
