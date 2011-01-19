cmake_minimum_required(VERSION 2.6)

project(TensorCalculus)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

add_definitions(-Wall)

include_directories(include)
include_directories(include/CC)
# include_directories(include/CG)
include_directories(include/DKTS)
include_directories(include/Matrix)
# include_directories(include/Mesh)
# include_directories(include/Quadrature)
include_directories(include/Representation)
# include_directories(include/Tensor)
# include_directories(include/Utilities)
# include_directories(include/Vector)

find_package(MKL)
if(MKL_FOUND)
  message(STATUS "Found MKL")
  set(BLAS_LIBRARIES ${MKL_LIBRARIES})
  set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
else(MKL_FOUND)
  # If MKL can't be found we're using the standard BLAS and LAPACK routines
  message(STATUS "Couldn't find MKL")
  find_package(BLAS REQUIRED)
  message(STATUS "Found BLAS")
  find_package(LAPACK REQUIRED)
  message(STATUS "Found LAPACK")
endif(MKL_FOUND)

#set(CMAKE_CXX_FLAGS "${CCMAKE_CXX_FLAGS} -fopenmp")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")



enable_testing()
add_subdirectory(test bin)
add_subdirectory(source)
