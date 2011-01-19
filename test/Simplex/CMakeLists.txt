add_executable(ValuesSimplexTest ValuesSimplexTest.cpp)
target_link_libraries(ValuesSimplexTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(SimplexTest SimplexTest.cpp)
target_link_libraries(SimplexTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(MeshTest MeshTest.cpp)
target_link_libraries(MeshTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(MeshFunctionTest MeshFunctionTest.cpp)
target_link_libraries(MeshFunctionTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(DomainTest DomainTest.cpp)
target_link_libraries(DomainTest ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
