#include <integrator_euler.hpp>
#include <integrator_RK2.hpp>
#include <iostream>
#include <math.h>

utils::Grid<double> my_rhs(utils::Grid<double> &waves,utils::Grid<double> &particles) {
	utils::Grid<double> rhs(waves);
	rhs += particles;
    return rhs;
}

utils::Grid<double> f_rhs_particles(utils::Grid<double> &waves,utils::Grid<double> &particles) {
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

utils::Grid<double> f_rhs_waves(utils::Grid<double> &waves,utils::Grid<double> &particles) {
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

//			if(ix==20 && iz==50) {
//				std::cout << " rhs " <<  rhs_local << "\n";
//			}

			// Add derivative of waves
			rhs_local += (waves.get_value(ix+1, iz) - 2.*waves.get_value(ix, iz)
					+ waves.get_value(ix-1, iz))/(del_x*del_x);
			rhs_local += (waves.get_value(ix, iz+1) - 2.*waves.get_value(ix, iz)
					+ waves.get_value(ix, iz-1))/(del_z*del_z);

			rhs_w.set_value(ix, iz, rhs_local);
//			if(ix==20 && iz==50) {
//				std::cout << " rhs " <<  rhs_w.get_value(ix, iz) << "\n";
//			}
		}
	}

//	std::cout << " end of rhs " << rhs_w.get_value(20, 50) << "\n";

	return rhs_w;
}


int main() {
	utils::Grid<double> waves(100,100);
	utils::Grid<double> particles(100,100);


	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_waves = f_rhs_waves;
	std::function<utils::Grid<double>(utils::Grid<double> &,utils::Grid<double> &)> rhs_particles = f_rhs_particles;


//	integrators::integrator_euler euler(10);
	integrators::integrator_RK2 euler(10);
	double time = 0.;
	double del_t = 0.001;
	del_t = 2.e-5;

	for(int i_time=0; i_time<50000; ++i_time) {
		euler.step(waves, particles, rhs_waves, rhs_particles, time, del_t);
		if(i_time%100==0) {
			utils::Grid<double> rhs = f_rhs_waves(waves, particles);
			std::cout << " step " << i_time << " ";
			std::cout << waves.get_value(10, 50) << " ";
			std::cout << waves.get_value(20, 50) << " ";
			std::cout << waves.get_value(30, 50) << " ";
			std::cout << rhs.get_value(30,50) << " ";
			std::cout << "\n";
			std::cout << "   --> ";
			std::cout << get_ana(waves,  10,  50) << " (";
			std::cout << get_ana(waves,  10,  50)-waves.get_value(10, 50) << "), ";
			std::cout << get_ana(waves,  20,  50) << " (";
			std::cout << get_ana(waves,  20,  50)-waves.get_value(20, 50) << "), ";
			std::cout << get_ana(waves,  30,  50) << " (";
			std::cout << get_ana(waves,  30,  50)-waves.get_value(30, 50) << ") ";
			std::cout << "\n";
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
