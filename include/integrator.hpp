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
  virtual ~integrator_base(){};

  virtual double step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
        double time, double del_t) = 0;
};

}  // namespace integrators

#endif
