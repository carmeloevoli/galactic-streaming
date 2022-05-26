#include "integrator_euler.hpp"

#include <iostream>

using namespace integrators;

integrator_euler::integrator_euler() {}

double integrator_euler::step(std::vector<utils::Grid<double>> &data,
    std::function<utils::Grid<double>(utils::Grid<double> &func_val)> rhs,
    double time, double del_t) {
  utils::Grid<double> rhs_grid = rhs(data[0]);
  std::cout << " Will do an Euler step - at some time in the near future...\n";
  time += del_t;
  return time;
}