add_executable(C4AmplitudenTest C4AmplitudenTest.cpp)
target_link_libraries(C4AmplitudenTest TensorCalculus)
target_link_libraries(C4AmplitudenTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(C4AmplitudenTest CC)
  
  
add_executable(C4Checks C4Checks.cpp)
target_link_libraries(C4Checks TensorCalculus)
target_link_libraries(C4Checks ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(C4Checks CC)
  
  
add_executable(C4Decomposition2 C4Decomposition2.cpp)
target_link_libraries(C4Decomposition2 TensorCalculus)
target_link_libraries(C4Decomposition2 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(C4Decomposition2 CC)
  
  
add_executable(C4MP2canonicalTest C4MP2canonicalTest.cpp)
target_link_libraries(C4MP2canonicalTest TensorCalculus)
target_link_libraries(C4MP2canonicalTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(C4MP2canonicalTest CC)
  
  
add_executable(C4MP2test C4MP2test.cpp)
target_link_libraries(C4MP2test TensorCalculus)
target_link_libraries(C4MP2test ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(C4MP2test CC)
  
  
add_executable(C4TransformationTest2 C4TransformationTest2.cpp)
target_link_libraries(C4TransformationTest2 TensorCalculus)
target_link_libraries(C4TransformationTest2 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(C4TransformationTest2 CC)


add_executable(MPSTest MPSTest.cpp)
target_link_libraries(MPSTest TensorCalculus)
target_link_libraries(MPSTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(MPSTest CC)
