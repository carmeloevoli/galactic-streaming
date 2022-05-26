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

  virtual double step(std::vector<utils::Grid<double>> &data, 
        std::function<utils::Grid<double>(utils::Grid<double> &)> ,
        double time, double del_t) = 0;
};

}  // namespace integrators

#endif