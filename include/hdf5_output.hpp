#ifndef HDF5_OUTPUT
#define HDF5_OUTPUT

#include "grid.hpp"

#include <string>
#include <hdf5.h>


namespace hdf5_IO {

class hdf5_writer{
public:
	hdf5_writer(std::string filename);
	int write_timestep(utils::Grid<double> &, utils::Grid<double> &, double d_time, int i_step);
	int write_grids(std::vector<double> x_grid, std::vector<double> z_grid);
private:
	hid_t open_file();
	void close_file();
	hid_t h5file, h5group;
	std::string m_filename;
	size_t i_entry;
};

}


#endif
