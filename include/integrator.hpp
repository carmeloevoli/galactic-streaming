#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <vector>

#include "grid.hpp"

/**
 * Generic integrator base class
 */
namespace integrators {

class integrator_base {
 public:
  virtual ~integrator_base(){};

  virtual double step(std::vector<utils::Grid<double>> &data, double time, double del_t) = 0;
};

}  // namespace integrators

#endif