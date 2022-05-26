cmake_minimum_required(VERSION 3.10)

add_executable (runit
	     main_hello.cpp  
  )

target_link_libraries (runit streaming ${STREAMING_EXTRA_LIBRARIES})

