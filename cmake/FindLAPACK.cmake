message(STATUS "In my own FindLAPACK.cmake")

find_library(LAPACK_LIBRARIES lapack PATHS ${LIB_INSTALL_DIR})
  
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARIES)

mark_as_advanced(LAPACK_LIBRARIES)


#message(STATUS "In my own FindLAPACK.cmake")
#
#find_library(LAPACK lapack PATHS ${LIB_INSTALL_DIR})
#  
#set(LAPACK_FOUND FALSE)
#
#if(LAPACK)
#  set(LAPACK_FOUND TRUE)
#  set(LAPACK_LIBRARIES ${LAPACK})
#  message(STATUE "Found LAPACK: ${LAPACK_LIBRARIES}")
#endif(LAPACK)
#
#mark_as_advanced(LAPACK_FOUND)
#mark_as_advanced(LAPACK_LIBRARIES)
