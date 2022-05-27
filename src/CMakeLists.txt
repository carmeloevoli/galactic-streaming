cmake_minimum_required(VERSION 3.10)

add_library(streaming SHARED
			"${git_revision_cpp}"
			hello_code.cpp
			integrator_euler.cpp
			params.cpp
			)


target_link_libraries(streaming ${STREAMING_EXTRA_LIBRARIES})

