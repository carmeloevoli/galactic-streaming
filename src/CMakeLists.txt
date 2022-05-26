cmake_minimum_required(VERSION 3.10)

add_library(streaming SHARED
    	hello_code.cpp
	"${git_revision_cpp}"
	)


target_link_libraries(streaming ${STREAMING_EXTRA_LIBRARIES})
