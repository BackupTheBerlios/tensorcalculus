add_executable(NewtonTest NewtonTest.cpp)
target_link_libraries(NewtonTest DKTS)
target_link_libraries(NewtonTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
