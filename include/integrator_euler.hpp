#ifndef INTEGRATOR_EULER_HPP
#define INTEGRATOR_EULER_HPP

#include "integrator.hpp"

namespace integrators {
class integrator_euler : public integrator_base {
 public:
  integrator_euler();

  double step(std::vector<utils::Grid<double>> &data,
        std::function<utils::Grid<double>(utils::Grid<double> &)> ,
        double time, double del_t);
};
}  // namespace integrators

#endif