cmake_minimum_required(VERSION 3.10)

add_library(streaming SHARED
	"${git_revision_cpp}"
	galaxy.cpp
	hello_code.cpp
	integrator.cpp
	integrator_euler.cpp
	integrator_RK2.cpp
	integrator_RK4.cpp
	shock_acceleration.cpp
	params.cpp
	utils.cpp
	)


target_link_libraries(streaming ${STREAMING_EXTRA_LIBRARIES})

