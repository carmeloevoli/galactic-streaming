#include "integrator_euler.hpp"

#include <iostream>

using namespace integrators;

integrator_euler::integrator_euler() {}

double integrator_euler::step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_waves,
		std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> &rhs_particles,
		double time, double del_t) {
	utils::Grid<double> rhs_grid = rhs_waves(waves,particles);
	std::cout << " Will do an Euler step - at some time in the near future...\n";
	time += del_t;
	return time;
}
