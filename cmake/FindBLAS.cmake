message(STATUS "In my own FindBLAS.cmake")

find_library(BLAS_LIBRARIES blas PATHS ${LIB_INSTALL_DIR})
  
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLAS DEFAULT_MSG BLAS_LIBRARIES)

mark_as_advanced(BLAS_LIBRARIES)
