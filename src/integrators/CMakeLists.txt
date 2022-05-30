cmake_minimum_required(VERSION 3.10)

add_library(integrators SHARED
	"${git_revision_cpp}"
	integrator.cpp
	integrator_euler.cpp
	integrator_RK2.cpp
	integrator_RK4.cpp
	)


target_link_libraries(integrators
 	PUBLIC ${STREAMING_EXTRA_LIBRARIES}
 )

