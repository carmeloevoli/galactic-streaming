#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "grid.hpp"

#include <vector>

/**
* Generic integrator base class
*/
namespace integrators{

class integrator_base {
    public:
    virtual ~integrator_base(){};

    virtual double step(std::vector<utils::Grid2D<double>> &data, double time, double del_t) = 0; 
};

}

#endif