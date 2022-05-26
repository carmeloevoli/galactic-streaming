#include <integrator_euler.hpp>
#include <iostream>

utils::Grid<double> my_rhs(utils::Grid<double> &grid) {
    return grid;
}


int main() {
  std::vector<utils::Grid<double>> grids;
  grids.push_back(utils::Grid<double>(100, 100));

std::function<utils::Grid<double>(utils::Grid<double> &)> func_rhs = my_rhs;


  integrators::integrator_euler euler;
  double time = 0.;
  double del_t = 0.1;
  euler.step(grids, my_rhs, time, del_t);

  return 0;
}