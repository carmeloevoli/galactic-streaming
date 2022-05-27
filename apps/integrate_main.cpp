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


	integrators::integrator_euler euler(10);
	double time = 0.;
	double del_t = 0.1;
	euler.step(waves, particles, rhs_waves, rhs_particles, time, del_t);

	std::cout << " Testing operators \n";
	utils::Grid<double> thingy(2,2,10.);
	std::cout << " At 1,1 " << thingy.get_value(1,1) << "\n";
	utils::Grid<double> other_thingy(2,2,1.);
	other_thingy = thingy;
	std::cout << " At 1,1 " << other_thingy.get_value(1,1) << "\n";
	other_thingy *= 2.4;
	std::cout << " At 1,1 " << other_thingy.get_value(1,1) << "\n";
	std::cout << " Before sum " << thingy.get_value(1,1) << " " << other_thingy.get_value(1, 1) << "\n";

	double val1=4., val2=2.;
	other_thingy = thingy*val1 + other_thingy*val2;
	std::cout << " should be 4*[0] + 2*[1] " << other_thingy.get_value(1, 1) << "\n";
	other_thingy += thingy;
	std::cout << other_thingy.get_value(1, 1) << " " << thingy.get_value(1,1) << "\n";

	return 0;
}
