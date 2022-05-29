#ifndef SHOCK_ACCELERATION_HPP
#define SHOCK_ACCELERATION_HPP

#include "grid.hpp"

#include <vector>


namespace shocks {
	class shock_acceleration{
public:
		shock_acceleration(int Nx, int Np);
		void make_coord_axes(int Nx, int Np);
		utils::Grid<double> get_rhs(utils::Grid<double> &waves, utils::Grid<double> &particles);
		utils::Grid<double> get_rhs_waves(utils::Grid<double> &waves, utils::Grid<double> &particles);
		std::vector<double> get_spatial_grid() const {return x_val;}
		std::vector<double> get_momentum_grid() const {return p_val;}
private:
		void set_velocity();
		std::vector<double> x_val, p_val;
		std::vector<double> velocity;
		std::vector<double> div_vel;
		double compression_ratio, del_x;
		double diffusion_strength, source_strength;
		double source_position = 0.;
		size_t m_Nx, m_Np;
		size_t index_sources;
	};
}

#endif
