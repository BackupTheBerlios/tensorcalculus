add_executable(LapackTest LapackTest.cpp)
target_link_libraries(LapackTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
