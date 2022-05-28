cmake_minimum_required(VERSION 3.10)

add_executable (runit
	     main_hello.cpp  
  )

target_link_libraries (runit streaming ${STREAMING_EXTRA_LIBRARIES})



add_executable (integrate
	     integrate_main.cpp  
  )

target_link_libraries (integrate streaming ${STREAMING_EXTRA_LIBRARIES})



add_executable (galaxy
	     galaxy_builder.cpp  
  )

target_link_libraries (galaxy streaming ${STREAMING_EXTRA_LIBRARIES})

add_executable (test_adaptive
	     test_integrate_adaptive.cpp  
  )

target_link_libraries (test_adaptive streaming ${STREAMING_EXTRA_LIBRARIES})
