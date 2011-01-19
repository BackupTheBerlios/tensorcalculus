add_executable(ACATest ACATest.cpp)
target_link_libraries(ACATest TensorCalculus)
target_link_libraries(ACATest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})              


add_executable(GramianMatrixInversionTest GramianMatrixInversionTest.cpp)
target_link_libraries(GramianMatrixInversionTest TensorCalculus)
target_link_libraries(GramianMatrixInversionTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

