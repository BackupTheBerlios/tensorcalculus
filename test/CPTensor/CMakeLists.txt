add_executable(CPTensorOpsTest CPTensorOpsTest.cpp)
target_link_libraries(CPTensorOpsTest ${BLAS_LIBRARIES})

add_executable(CPTensorTest CPTensorTest.cpp)
add_executable(TensorBasisOpsTest TensorBasisOpsTest.cpp)
add_executable(TensorBasisTest TensorBasisTest.cpp)
