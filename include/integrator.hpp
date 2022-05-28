#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <vector>
#include <functional>

#include "grid.hpp"

/**
 * Generic integrator base class
 */
namespace integrators {

class integrator_base {
 public:
	integrator_base(int verbosity) {
	  m_verbosity = verbosity;
	  scheme_order = 1.;
	  beta = 0.8; // Safety factor for adaptive step-size control
	}
  virtual ~integrator_base(){
  };

  virtual double step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
        double time, double del_t) = 0;
  double step_adaptive(utils::Grid<double> &waves,  utils::Grid<double> &particles,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
        double time, double & del_t, double epsilon);
 protected:
  int m_verbosity;
  double scheme_order, beta;
};


}  // namespace integrators

#endif
