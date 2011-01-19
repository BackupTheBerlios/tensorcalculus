add_executable(VectorOpsTest VectorOpsTest.cpp)
target_link_libraries(VectorOpsTest ${BLAS_LIBRARIES})

add_executable(WeakVectorOpsTest WeakVectorOpsTest.cpp)
target_link_libraries(WeakVectorOpsTest ${BLAS_LIBRARIES})

add_executable(WeakVectorTest WeakVectorTest.cpp)
target_link_libraries(WeakVectorTest ${BLAS_LIBRARIES})
