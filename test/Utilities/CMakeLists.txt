add_executable(FactorialsTest FactorialsTest.cpp)

add_executable(IndexTest IndexTest.cpp)
target_link_libraries(IndexTest TensorCalculus)

add_executable(SkippedProductsTest SkippedProductsTest.cpp)
target_link_libraries(SkippedProductsTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
