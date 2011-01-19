# Find all source files in the project tree
# Attention: This isn't optimal because files which aren't
#            officially part of the project will be added too!
file(GLOB_RECURSE ALL_SOURCE_FILES ${CMAKE_SOURCE_DIR}/source/CC/*.cpp ${CMAKE_SOURCE_DIR}/source/CC/*.c)

# Find all header files in the include directory
file(GLOB ALL_HEADER_FILES ${CMAKE_SOURCE_DIR}/include/CC/*.hpp ${CMAKE_SOURCE_DIR}/include/CC/*.h)

# Add the CC library
# Adding ALL_HEADER_FILES to TensorCalculus is nescessary due to the MSVS project file
add_library(CC ${ALL_SOURCE_FILES} ${ALL_HEADER_FILES})

target_link_libraries(CC TensorCalculus)
