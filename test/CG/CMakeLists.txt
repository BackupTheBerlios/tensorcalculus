add_executable(CGMethodTest CGMethodTest.cpp)
target_link_libraries(CGMethodTest ${BLAS_LIBRARIES})

add_executable(ThreePGTest ThreePGTest.cpp)
add_executable(BisectionTest BisectionTest.cpp)

add_executable(NetworkCGMethodTest NetworkCGMethodTest.cpp)
target_link_libraries(NetworkCGMethodTest ${BLAS_LIBRARIES})
target_link_libraries(NetworkCGMethodTest TensorCalculus)
