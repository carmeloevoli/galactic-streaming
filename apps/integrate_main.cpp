#include <integrator_euler.hpp>
#include <iostream>

utils::Grid<double> my_rhs(utils::Grid<double> &waves,utils::Grid<double> &particles) {
	utils::Grid<double> rhs(waves);
	rhs += particles;
    return rhs;
}


int main() {
	utils::Grid<double> waves(100,100);
	utils::Grid<double> particles(100,100);


	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_waves = my_rhs;
	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_particles = my_rhs;


	integrators::integrator_euler euler;
	double time = 0.;
	double del_t = 0.1;
	euler.step(waves, particles, rhs_waves, rhs_particles, time, del_t);

	return 0;
}
