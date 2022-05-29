#include "shock_acceleration.hpp"

#include "utils.hpp"

#include <iostream>

namespace shocks {

shock_acceleration::shock_acceleration(int Nx, int Np) {

	compression_ratio = 4.;
	make_coord_axes(Nx, Np);
	diffusion_strength = 0.1;
	source_strength = 1.;

}

void shock_acceleration::make_coord_axes(int Nx, int Np) {
	m_Nx = Nx;
	m_Np = Np;
	// Make linear grid in space (normalised coordinate from 0 to 1)
	x_val = utils::build_lin_axis<double>(-0.5, 0.5, Nx);
	p_val = utils::build_log_axis<double>(1., 100., Np);

	del_x = x_val[1] - x_val[0];

	auto it = find(x_val.begin(), x_val.end(), source_position);

	if (it != x_val.end())
	{

		// calculating the index
		// of K
		index_sources = it - x_val.begin();
		std::cout << " Source located at " << index_sources << "\n";
	} else {
		index_sources = -1;
		std::cerr<<" Source position not on spatial grid \n";
	}

	set_velocity();
}

void shock_acceleration::set_velocity() {
	velocity.resize(m_Nx);

	 for(size_t i_x=0; i_x < velocity.size(); ++i_x) {
		 if(x_val[i_x] < 0.) {
			 velocity[i_x] = 1.;
		 } else {
			 velocity[i_x] = 1/compression_ratio;
		 }
	 }

	 // set divergence to zero everwhere
	 div_vel.resize(m_Nx);
	 std::fill(div_vel.begin(), div_vel.end(), 0);

	 // Compute divergence of velocity using finite differences
	 for(size_t i_x=1; i_x < velocity.size()-1; ++i_x) {
		 div_vel[i_x] = (velocity[i_x+1] - velocity[i_x-1])/(x_val[i_x+1] - x_val[i_x-1]);
	 }

}

utils::Grid<double> shock_acceleration::get_rhs(utils::Grid<double> &waves,
		utils::Grid<double> &particles) {
	// Waves are not needed in this test case
	utils::Grid<double> rhs(particles);
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();

	assert(Nx==m_Nx && Nz==m_Np);

	// Loop over both grids
	for(size_t i_x=1; i_x<Nx-1; ++i_x) {
		for(size_t i_p=1; i_p<Nz-1; ++i_p) {
			double rhs_local = 0.;

			// Contribution of source term
			if(i_x == index_sources && i_p==2) {
				rhs_local += source_strength;
			}

			// Contribution of spatial evolution
			rhs_local -= particles.get(i_x, i_p)*div_vel[i_x];
			double grad_j = (particles.get(i_x+1,i_p) - particles.get(i_x-1,i_p))/(2.*del_x);
			rhs_local -= velocity[i_x]*(grad_j);
			double laplace_j = (particles.get(i_x+1,i_p) - 2*particles.get(i_x,i_p)
					+ particles.get(i_x-1,i_p))/(del_x*del_x);

			rhs_local += diffusion_strength*laplace_j;

			// Contribution of momentum evolution
			double mom_deriv_pj = (particles.get(i_x,i_p+1)*p_val[i_p+1] -
					particles.get(i_x,i_p-1)*p_val[i_p-1])/(p_val[i_p+1] - p_val[i_p-1]);
			rhs_local += 1./3.*div_vel[i_x]*mom_deriv_pj;

			rhs.set_value(i_x, i_p, rhs_local);
		}
	}

	return rhs;

}

utils::Grid<double> shock_acceleration::get_rhs_waves(utils::Grid<double> &waves, utils::Grid<double> &particles) {
	// Waves are not needed in this test case
	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	utils::Grid<double> rhs(Nx, Nz, 0.);
	return rhs;
}


}
