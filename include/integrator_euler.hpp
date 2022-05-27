#ifndef INTEGRATOR_EULER_HPP
#define INTEGRATOR_EULER_HPP

#include "integrator.hpp"

namespace integrators {
class integrator_euler : public integrator_base {
public:
	integrator_euler(int verbosity);

	double step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
			std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
			std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
			double time, double del_t);
};
}  // namespace integrators

#endif
