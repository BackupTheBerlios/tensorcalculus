add_executable(Eigenfunctions Eigenfunctions.cpp)
target_link_libraries(Eigenfunctions ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
