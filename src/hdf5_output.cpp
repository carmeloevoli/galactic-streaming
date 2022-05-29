#include "hdf5_output.hpp"

#include <sstream>

using namespace hdf5_IO;


hdf5_writer::hdf5_writer(std::string filename) {
	m_filename = filename;
	h5file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	h5group  =  H5Gcreate2(h5file, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	i_entry = 0;

	close_file();
}


hid_t hdf5_writer::open_file() {
	h5file = H5Fopen(m_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	h5group = H5Gopen2(h5file, "Data", H5P_DEFAULT);
	return h5group;
}


void hdf5_writer::close_file() {
	H5Gclose(h5group);
	H5Fclose(h5file);
}

int hdf5_writer::write_grids(std::vector<double> x_grid, std::vector<double> z_grid) {

	open_file();

	hid_t datatype_grid = H5Tcopy(H5T_NATIVE_DOUBLE);

	// Create the dataspace
	hsize_t DimsData = x_grid.size();
	hid_t dataspace_x = H5Screate_simple(1, &DimsData, NULL);

	// Save position data
	//	std::cout << " Trying to create " << name_pos.c_str() << std::endl;
	hid_t dataset_x_grid = H5Dcreate2(h5group, "x_pos", datatype_grid, dataspace_x,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_x_grid, datatype_grid, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_grid.data());
	H5Dclose(dataset_x_grid);

	// Create the dataspace
	DimsData = z_grid.size();
	hid_t dataspace_z = H5Screate_simple(1, &DimsData, NULL);

	// Save position data
	//	std::cout << " Trying to create " << name_pos.c_str() << std::endl;
	hid_t dataset_z_grid = H5Dcreate2(h5group, "z_pos", datatype_grid, dataspace_z,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_z_grid, datatype_grid, H5S_ALL, H5S_ALL, H5P_DEFAULT, z_grid.data());
	H5Dclose(dataset_z_grid);

	close_file();

	return 0;
}


int hdf5_writer::write_timestep(utils::Grid<double> &waves, utils::Grid<double> &particles,
	double d_time, int i_step) {
	// Reopen file
	open_file();

	std::ostringstream oss;
	oss << "step" << i_step;
	std::string gname = oss.str();
	hid_t h5local_group  =  H5Gcreate2(h5group, gname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create datatype
	hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);


	// add time to local group
	hid_t AttrSpace = H5Screate(H5S_SCALAR);
	hid_t attr_time = H5Acreate2(h5local_group, "time", datatype, AttrSpace,
				H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_time, datatype, &d_time);

	size_t Nx=waves.get_Nx();
	size_t Nz=waves.get_Nz();
	hsize_t Nx_grid[2] = {Nx, Nz};

	// Create new dataspace
	hid_t dataspace2D = H5Screate_simple(2, Nx_grid, NULL);


	// Save waves data
	hid_t dataset_flux = H5Dcreate2(h5local_group, "waves", datatype, dataspace2D,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	std::vector<double> raw_data = waves.get_data_copy();
	H5Dwrite(dataset_flux, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw_data.data());
	H5Dclose(dataset_flux);

	dataset_flux = H5Dcreate2(h5local_group, "particles", datatype, dataspace2D,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	raw_data = particles.get_data_copy();
	H5Dwrite(dataset_flux, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw_data.data());
	H5Dclose(dataset_flux);


	// Close dataspace
	H5Sclose(dataspace2D);

	i_entry++;

	H5Gclose(h5local_group);

	// Close file again
	close_file();


	return 0;
}
