#include "integrator_euler.hpp"

#include <iostream>

using namespace integrators;

integrator_euler::integrator_euler() {}

double integrator_euler::step(std::vector<utils::Grid<double>> &data, double time, double del_t) {
  std::cout << " Will do an Euler step - at some time in the near future...\n";
  time += del_t;
  return time;
}