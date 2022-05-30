#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "grid.hpp"
#include "coupled_pde.hpp"

#include <vector>
#include <functional>
#include <memory>


/**
 * Generic integrator base class
 */
namespace integrators {

class integrator_base {
 public:
	integrator_base(std::shared_ptr<coupled_pdes::coupled_pde> pde, int verbosity) : m_pde(pde) {
		m_verbosity = verbosity;
		scheme_order = 1.;
		beta = 0.8; // Safety factor for adaptive step-size control
	}
	integrator_base(int verbosity) {
		m_verbosity = verbosity;
		scheme_order = 1.;
		beta = 0.8; // Safety factor for adaptive step-size control
	}
	virtual ~integrator_base(){
	};

//  virtual double step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
//        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
//        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
//        double time, double del_t) = 0;
  virtual double step(utils::Grid<double> &waves,  utils::Grid<double> &particles,
        double time, double del_t) = 0;
//  double step_adaptive(utils::Grid<double> &waves,  utils::Grid<double> &particles,
//        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
//        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
//        double time, double & del_t, double epsilon);
  double step_adaptive(utils::Grid<double> &waves,  utils::Grid<double> &particles,
        double time, double & del_t, double epsilon);
 protected:
  std::shared_ptr<coupled_pdes::coupled_pde> m_pde;
  int m_verbosity;
  double scheme_order, beta;
};


}  // namespace integrators

#endif
