#include <iostream>
#include <integrator_euler.hpp>

int main() {

    std::vector<utils::Grid2D<double>> grids;
    grids.push_back(utils::Grid2D<double>(100, 100));

    integrators::integrator_euler euler;
    double time = 0.;
    double del_t = 0.1;
    euler.step(grids, time, del_t);

    return 0;
}