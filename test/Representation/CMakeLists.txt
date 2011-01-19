add_executable(ComponentTensorTest ComponentTensorTest.cpp)
target_link_libraries(ComponentTensorTest TensorCalculus)
target_link_libraries(ComponentTensorTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(ComponentTensorTest DKTS)
target_link_libraries(ComponentTensorTest CC)


add_executable(CPALSTest CPALSTest.cpp)
target_link_libraries(CPALSTest TensorCalculus)
target_link_libraries(CPALSTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


add_executable(CPTensorRepresentationTest CPTensorRepresentationTest.cpp)
target_link_libraries(CPTensorRepresentationTest ${BLAS_LIBRARIES})
target_link_libraries(CPTensorRepresentationTest TensorCalculus)
target_link_libraries(CPTensorRepresentationTest ${LAPACK_LIBRARIES})
target_link_libraries(CPTensorRepresentationTest DKTS)
  

add_executable(TensorChainTest TensorChainTest.cpp)
target_link_libraries(TensorChainTest TensorCalculus)
target_link_libraries(TensorChainTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})


add_executable(TensorRepresentationHeaderTest TensorRepresentationHeaderTest.cpp)
target_link_libraries(TensorRepresentationHeaderTest ${BLAS_LIBRARIES})
target_link_libraries(TensorRepresentationHeaderTest TensorCalculus)

