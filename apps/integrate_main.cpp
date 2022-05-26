#include <integrator_euler.hpp>
#include <iostream>

int main() {
  std::vector<utils::Grid<double>> grids;
  grids.push_back(utils::Grid<double>(100, 100));

  integrators::integrator_euler euler;
  double time = 0.;
  double del_t = 0.1;
  euler.step(grids, time, del_t);

  return 0;
}