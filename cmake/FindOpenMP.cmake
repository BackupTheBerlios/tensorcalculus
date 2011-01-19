message(STATUS "In my own FindOpenMP.cmake")

find_library(GOMP_LIBRARIES gomp PATHS ${LIB_INSTALL_DIR})
  
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GOMP DEFAULT_MSG GOMP_LIBRARIES)

mark_as_advanced(GOMP_LIBRARIES)
