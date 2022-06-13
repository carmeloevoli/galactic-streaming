cmake_minimum_required(VERSION 3.1)
project(streaming 
	LANGUAGES C CXX)
set(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(STREAMING_EXTRA_SOURCES)
set(STREAMING_EXTRA_INCLUDES)
set(STREAMING_EXTRA_LIBRARIES)

set(CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

# specify the C standard
set(C_STANDARD 14)
set(CMAKE_C_STANDARD_REQUIRED True)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -pedantic -ftree-vectorize -O3")
set(STREAMING_INCLUDES)
set(STREAMING_LIBRARIES)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# PLOG (provided)
include_directories(external/plog/include)

# C++ Threads required
#find_package(Threads REQUIRED)
#list(APPEND STREAMING_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# GSL (required)
find_package(GSL REQUIRED)
list(APPEND STREAMING_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND STREAMING_EXTRA_LIBRARIES ${GSL_LIBRARIES})

# Version info from Git
option(ENABLE_GIT "Embedding information about code version from git" ON)
if(ENABLE_GIT)
    include(GetGitRevisionDescription)
    get_git_head_revision(GIT_REFSPEC GIT_SHA1)
    git_local_changes(GIT_HAS_LOCAL_CHANGES)
    set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
    configure_file("src/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)
endif(ENABLE_GIT)
IF(NOT ENABLE_GIT OR (${GIT_SHA1} STREQUAL "GITDIR-NOTFOUND"))
    set(GIT_REFSPEC "")
    set(GIT_SHA1 "")
    set(GIT_DESC "${STREAMING_RELEASE_VERSION}-no-git")
endif()
message(STATUS "Code version: ${GIT_DESC} ${GIT_SHA1} ${GIT_REFSPEC}")

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
include_directories(include ${STREAMING_EXTRA_INCLUDES})

#add_library(streaming SHARED
#    src/main.cpp
#    "${git_revision_cpp}"
#    )
#    target_link_libraries(streaming ${STREAMING_EXTRA_LIBRARIES})

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# copy data files in build
#file(GLOB STREAMING_DATA "data/*.txt")
#file(COPY ${STREAMING_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )

find_package(HDF5 REQUIRED)
MESSAGE(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIR}")
# Make sure serial hdf5 library works in serial form
add_library(SerialHdf5 INTERFACE)
target_link_libraries(SerialHdf5
  INTERFACE ${HDF5_LIBRARIES}
  )

if(${HDF5_IS_PARALLEL})
  
  find_package(MPI)
  MESSAGE(" Loaded MPI for parallel HDF5 ")
  MESSAGE(STATUS " Found include path at ${MPI_CXX_INCLUDE_PATH}")

  target_include_directories(SerialHdf5
    INTERFACE ${HDF5_INCLUDE_DIR}
    INTERFACE ${MPI_CXX_INCLUDE_PATH}
    )
  target_compile_definitions(SerialHdf5 INTERFACE OMPI_SKIP_MPICXX)

else()

  target_include_directories(SerialHdf5
    INTERFACE ${HDF5_INCLUDE_DIR}
    )
  
endif()


# make examples
add_subdirectory(src)
add_subdirectory(apps)

#add_executable(runit src/main.cpp)

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "STREAMING_INCLUDES ...... = ${STREAMING_EXTRA_INCLUDES}")
message (STATUS "STREAMING_LIBRARIES ..... = ${STREAMING_EXTRA_LIBRARIES}")
