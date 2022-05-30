cmake_minimum_required(VERSION 3.10)

# add integrators
add_subdirectory(integrators)


add_library(streaming SHARED
	"${git_revision_cpp}"
	galaxy.cpp
	hello_code.cpp
	shock_acceleration.cpp
	params.cpp
	utils.cpp
	hdf5_output.cpp
	coupled_pde.cpp
	)


target_link_libraries(streaming
	PUBLIC integrators
 	PUBLIC ${STREAMING_EXTRA_LIBRARIES}
 	PUBLIC SerialHdf5
 )

