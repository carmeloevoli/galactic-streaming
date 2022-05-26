#ifndef INTEGRATOR_EULER_HPP
#define INTEGRATOR_EULER_HPP

#include "integrator.hpp"

namespace integrators {
    class integrator_euler : public integrator_base {
    public:
        integrator_euler();

        double step(std::vector<utils::Grid2D<double>> &data, double time, double del_t); 


    };
}

#endif