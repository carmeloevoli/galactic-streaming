#include "integrators/integrator_euler.hpp"
#include "integrators/integrator_RK2.hpp"
#include "integrators/integrator_RK4.hpp"

#include <iostream>
#include <math.h>


class test_time_dependence : public coupled_pdes::coupled_pde {
public:
	test_time_dependence(){}

	utils::Grid<double> get_rhs_waves(utils::Grid<double> &waves,
			utils::Grid<double> &particles) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();

		std::vector<double> xGrid(Nx);
		std::vector<double> zGrid(Nz);

		double del_x = 1./(Nx-1.);
		double del_z = 1./(Nz-1.);
		for(size_t ix=0; ix<Nx; ++ix) {
			xGrid[ix] = ix*del_x;
		}
		for(size_t iz=0; iz<Nz; ++iz) {
			zGrid[iz] = iz*del_z;
		}
		utils::Grid<double> rhs_w(waves);
		rhs_w.clear_grid();
		for(size_t ix=1; ix<Nx-1; ++ix) {
			for(size_t iz=1; iz<Nz-1; ++iz) {
				// Add derivative of waves
				double rhs_local = (waves.get_value(ix+1, iz) - 2.*waves.get_value(ix, iz)
						+ waves.get_value(ix-1, iz))/(del_x*del_x);
				rhs_local += (waves.get_value(ix, iz+1) - 2.*waves.get_value(ix, iz)
						+ waves.get_value(ix, iz-1))/(del_z*del_z);

				rhs_w.set_value(ix, iz, rhs_local);
			}
		}

		return rhs_w;
	}

	utils::Grid<double> get_rhs_particles(utils::Grid<double> &waves,
			utils::Grid<double> &particles) {
		utils::Grid<double> rhs_p(particles);
		return rhs_p;
	}

	double get_ana_time_dep(utils::Grid<double> &waves, double time, int t_ix, int t_iz) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();

		std::vector<double> xGrid(Nx);
		std::vector<double> zGrid(Nz);

		double del_x = 1./(Nx-1.);
		double del_z = 1./(Nz-1.);
		for(size_t ix=0; ix<Nx; ++ix) {
			xGrid[ix] = ix*del_x;
		}
		for(size_t iz=0; iz<Nz; ++iz) {
			zGrid[iz] = iz*del_z;
		}

		double xVal = xGrid[t_ix];
		double zVal = zGrid[t_iz];
		double pi = M_PI;

		return sin(pi*xVal)*sin(pi*zVal)*exp(-2*pi*time);

	}

	void set_ana_time_dep(utils::Grid<double> &waves, double time) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();
		for(size_t ix=1; ix<Nx-1; ++ix) {
			for(size_t iz=1; iz<Nz-1; ++iz) {
				waves.set_value(ix, iz, get_ana_time_dep(waves, time, ix, iz));
			}
		}
	}

	double get_l2error_time(utils::Grid<double> &waves, double time) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();

		double sum_diff_sqr=0;
		double norm=0;
		for(size_t ix=1; ix<Nx-1; ++ix) {
			for(size_t iz=1; iz<Nz-1; ++iz) {
				double ana = get_ana_time_dep(waves, time, ix, iz);
				double num = waves.get_value(ix, iz);
				double diffsqr = (ana-num)*(ana-num);
				sum_diff_sqr += diffsqr;
				norm += ana*ana;
			}
		}
		double l2err = sqrt(sum_diff_sqr/norm);
		return l2err;
	}

};

class test_time_exp: public coupled_pdes::coupled_pde {
public:
	test_time_exp(){}

	utils::Grid<double> get_rhs_waves(utils::Grid<double> &waves,utils::Grid<double> &particles) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();
		double lambd = -2.;
		utils::Grid<double> rhs_w(Nx, Nz);
		rhs_w = waves*lambd;
		return rhs_w;
	}

	utils::Grid<double> get_rhs_particles(utils::Grid<double> &waves,
			utils::Grid<double> &particles) {
		utils::Grid<double> rhs_p(particles);
		return rhs_p;
	}

};

class test_steady_state : public coupled_pdes::coupled_pde {
public:

	test_steady_state() {}

	utils::Grid<double> get_rhs_waves(utils::Grid<double> &waves,utils::Grid<double> &particles) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();

		std::vector<double> xGrid(Nx);
		std::vector<double> zGrid(Nz);

		double del_x = 1./(Nx-1.);
		double del_z = 1./(Nz-1.);
		for(size_t ix=0; ix<Nx; ++ix) {
			xGrid[ix] = ix*del_x;
		}
		for(size_t iz=0; iz<Nz; ++iz) {
			zGrid[iz] = iz*del_z;
		}
		double pi = M_PI;
		utils::Grid<double> rhs_w(waves);
		rhs_w.clear_grid();
		for(size_t ix=1; ix<Nx-1; ++ix) {
			double xVal = xGrid[ix];
			for(size_t iz=1; iz<Nz-1; ++iz) {
				double zVal = zGrid[iz];

				double rhs_local =2.*pi*pi*sin(pi*xVal)*sin(pi*zVal);

				// Add derivative of waves
				rhs_local += (waves.get_value(ix+1, iz) - 2.*waves.get_value(ix, iz)
						+ waves.get_value(ix-1, iz))/(del_x*del_x);
				rhs_local += (waves.get_value(ix, iz+1) - 2.*waves.get_value(ix, iz)
						+ waves.get_value(ix, iz-1))/(del_z*del_z);

				rhs_w.set_value(ix, iz, rhs_local);
			}
		}

		return rhs_w;
	}

	utils::Grid<double> get_rhs_particles(utils::Grid<double> &waves,utils::Grid<double> &particles) {
		utils::Grid<double> rhs_p(particles);
		return rhs_p;
	}

	double get_ana(utils::Grid<double> &waves, int t_ix, int t_iz) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();

		std::vector<double> xGrid(Nx);
		std::vector<double> zGrid(Nz);

		double del_x = 1./(Nx-1.);
		double del_z = 1./(Nz-1.);
		for(size_t ix=0; ix<Nx; ++ix) {
			xGrid[ix] = ix*del_x;
		}
		for(size_t iz=0; iz<Nz; ++iz) {
			zGrid[iz] = iz*del_z;
		}

		double xVal = xGrid[t_ix];
		double zVal = zGrid[t_iz];
		double pi = M_PI;

		return sin(pi*xVal)*sin(pi*zVal);

	}

	double get_l2error(utils::Grid<double> &waves) {
		size_t Nx=waves.get_Nx();
		size_t Nz=waves.get_Nz();

		double sum_diff_sqr=0;
		double norm=0;
		for(size_t ix=1; ix<Nx-1; ++ix) {
			for(size_t iz=1; iz<Nz-1; ++iz) {
				double ana = get_ana(waves, ix, iz);
				double num = waves.get_value(ix, iz);
				double diffsqr = (ana-num)*(ana-num);
				sum_diff_sqr += diffsqr;
				norm += ana*ana;
			}
		}
		double l2err = sqrt(sum_diff_sqr/norm);
		return l2err;
	}



};

enum class TestType {
	SteadyState,
	TimeDependent,
	TimeExp
};

int main() {
	utils::Grid<double> waves(100,100);
	utils::Grid<double> particles(100,100);

	double time = 0.;
	double del_t = 0.001;
	int i_time = 0;
	del_t = 2.e-6;

	TestType my_test = TestType::TimeExp;

	if(my_test == TestType::SteadyState) {

		std::shared_ptr<test_steady_state> steady = std::make_shared<test_steady_state>();

		//	integrators::integrator_euler integrator(steady, 10);
		//	integrators::integrator_RK2 integrator(steady, 10);
		integrators::integrator_RK4 integrator(steady, waves,10);

		//	for(int i_time=0; i_time<50000; ++i_time) {
		while(time<1.) {
			time = integrator.step(waves, particles, time, del_t);
			if(i_time%100==0) {
				utils::Grid<double> rhs = steady->get_rhs_waves(waves, particles);
				std::cout << " step " << i_time << " -> " << time << " ";
				std::cout << waves.get_value(10, 50) << " ";
				std::cout << waves.get_value(20, 50) << " ";
				std::cout << waves.get_value(30, 50) << " ";
				std::cout << rhs.get_value(30,50) << " ";
				std::cout << "\n";
				std::cout << "   --> ";
				std::cout << steady->get_ana(waves,  10,  50) << " (";
				std::cout << steady->get_ana(waves,  10,  50)-waves.get_value(10, 50) << "), ";
				std::cout << steady->get_ana(waves,  20,  50) << " (";
				std::cout << steady->get_ana(waves,  20,  50)-waves.get_value(20, 50) << "), ";
				std::cout << steady->get_ana(waves,  30,  50) << " (";
				std::cout << steady->get_ana(waves,  30,  50)-waves.get_value(30, 50) << ") ";
				std::cout << "\n";
			}
			i_time++;
		}

		std::cout << " l2 error: " << steady->get_l2error(waves) << "\n";

	} else if (my_test== TestType::TimeDependent) {

		std::shared_ptr<test_time_dependence> time_dep = std::make_shared<test_time_dependence>();

		integrators::integrator_RK4 integrator(time_dep, waves,10);

		double time = 0.;
		double del_t = 0.001;
		int i_time = 0;
		del_t = 2.e-5;

		time_dep->set_ana_time_dep(waves, 0.);
		while(time<0.01) {

			time = integrator.step(waves, particles, time, del_t);
			if(i_time%100==0) {
				std::cout << " step " << i_time << " at " << time << " -> ";
				std::cout << waves.get_value(10, 50) << " ";
				std::cout << waves.get_value(20, 50) << " ";
				std::cout << waves.get_value(30, 50) << " ";
				std::cout << "\n";
				std::cout << "   --> ";
				std::cout << time_dep->get_ana_time_dep(waves, time, 10,  50) << " (";
				std::cout << time_dep->get_ana_time_dep(waves, time, 10,  50)-waves.get_value(10, 50) << "), ";
				std::cout << time_dep->get_ana_time_dep(waves, time, 20,  50) << " (";
				std::cout << time_dep->get_ana_time_dep(waves, time, 20,  50)-waves.get_value(20, 50) << "), ";
				std::cout << time_dep->get_ana_time_dep(waves, time, 30,  50) << " (";
				std::cout << time_dep->get_ana_time_dep(waves, time, 30,  50)-waves.get_value(30, 50) << ") ";
				std::cout << "\n";
				double l2_error = time_dep->get_l2error_time(waves, time);
				std::cout << " l2 error " << l2_error << "\n";
			}
			i_time++;

		}
	} else if(my_test==TestType::TimeExp) {
		double time = 0.;
		double del_t = 0.001;
		int i_time = 0;
		del_t = 0.01;

		std::shared_ptr<test_time_exp> time_exp = std::make_shared<test_time_exp>();

		integrators::integrator_RK4 integrator(time_exp, waves,10);


		waves.set_value(2, 2, 1.);

		while(time<1.) {

			time = integrator.step(waves, particles, time, del_t);
			if(i_time%10==0) {
				std::cout << " step " << i_time << " at " << time << " -> ";
				std::cout << waves.get_value(2, 2) << " ";
				std::cout << waves.get_value(2, 2) - exp(-2.*time);
				std::cout << "\n";
			}
			i_time++;

		}
	}

//	std::cout << " Testing operators \n";
//	utils::Grid<double> thingy(2,2,10.);
//	std::cout << " At 1,1 " << thingy.get_value(1,1) << "\n";
//	utils::Grid<double> other_thingy(2,2,1.);
//	other_thingy = thingy;
//	std::cout << " At 1,1 " << other_thingy.get_value(1,1) << "\n";
//	other_thingy *= 2.4;
//	std::cout << " At 1,1 " << other_thingy.get_value(1,1) << "\n";
//	std::cout << " Before sum " << thingy.get_value(1,1) << " " << other_thingy.get_value(1, 1) << "\n";
//
//	double val1=4., val2=2.;
//	other_thingy = thingy*val1 + other_thingy*val2;
//	std::cout << " should be 4*[0] + 2*[1] " << other_thingy.get_value(1, 1) << "\n";
//	other_thingy += thingy;
//	std::cout << other_thingy.get_value(1, 1) << " " << thingy.get_value(1,1) << "\n";
//
//	utils::Grid<double> thingy2(100,100,1.);
//	thingy2.set_value(20,50, 11.7021);
//	std::cout << " thingy2 value " << thingy2.get_value(20, 50) << " " << thingy2.get_size()<< "\n";
//	utils::Grid<double> test_thingy = thingy2*0.1;
//	std::cout << " test thingy  " << test_thingy.get_value(20, 50) << "\n";
//	thingy2 += thingy2*2.2;
////	std::cout << " final " << thingy2.get_value(20,50) << " " << thingy2.get_value(22,50) << "\n";
//
//	utils::Grid<double> rhs = f_rhs_waves(waves, particles);
//	std::cout << " rhs value " << rhs.get_value(20, 50) << " " << rhs.get_size() << "\n";
//	utils::Grid<double> test_me = rhs*0.1;
//	std::cout << " test value " << test_me.get_value(20, 50) << "\n";


	return 0;
}
