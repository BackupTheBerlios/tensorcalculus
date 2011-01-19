# Find all source files in the project tree
# Attention: This isn't optimal because files which aren't
#            officially part of the project will be added too!
file(GLOB_RECURSE ALL_SOURCE_FILES ${CMAKE_SOURCE_DIR}/source/DKTS/*.cpp ${CMAKE_SOURCE_DIR}/source/DKTS/*.c)

# Find all header files in the include directory
file(GLOB ALL_HEADER_FILES ${CMAKE_SOURCE_DIR}/include/DKTS/*.hpp ${CMAKE_SOURCE_DIR}/include/DKTS/*.h)

# Add the TensorCalculus library
# Adding ALL_HEADER_FILES to TensorCalculus is nescessary due to the MSVS project file
add_library(DKTS ${ALL_SOURCE_FILES} ${ALL_HEADER_FILES})

# Grouping of the source files and their associated headers in 
# the Visual Studio project file
foreach(SOURCE_FILE ${ALL_SOURCE_FILES})
  get_filename_component(HEADER_FILE ${SOURCE_FILE} NAME_WE)
  set(HEADER_FILE ${CMAKE_SOURCE_DIR}/include/DKTS/${HEADER_FILE}.hpp)
  list(REMOVE_ITEM ALL_HEADER_FILES ${HEADER_FILE})
  #message(STATUS ${HEADER_FILE})
  string(REGEX REPLACE ${CMAKE_CURRENT_SOURCE_DIR} "" REL_PATH "${SOURCE_FILE}")
  string(REGEX REPLACE "[\\\\/][^\\\\/]*$" "" REL_PATH "${REL_PATH}")
  string(REGEX REPLACE "^[\\\\/]" "" REL_PATH "${REL_PATH}")
  string(REGEX REPLACE "/" "\\\\\\\\" REL_PATH "${REL_PATH}")
  # if("${REL_PATH}" STREQUAL "")
    # message(STATUS "Source" " FILES: " ${SOURCE_FILE})
    # message(STATUS "Header" " FILES: " ${HEADER_FILE})
    # source_group("Source" FILES ${SOURCE_FILE})
    # source_group("Header" FILES ${HEADER_FILE})
  # else("${REL_PATH}" STREQUAL "")
    # message(STATUS "Source Files\\\\${REL_PATH} FILES ${SOURCE_FILE}")
    # message(STATUS "Header Files\\\\${REL_PATH} FILES ${HEADER_FILE}")
    source_group("Source Files\\\\${REL_PATH}" FILES ${SOURCE_FILE})
    source_group("Header Files\\\\${REL_PATH}" FILES ${HEADER_FILE})
  # endif("${REL_PATH}" STREQUAL "")
endforeach(SOURCE_FILE)

# message("Remaining header files:")

foreach(HEADER_FILE ${ALL_HEADER_FILES})
  # message(STATUS "Header Files\\\\Base FILES ${HEADER_FILE}")
  source_group("Header Files\\\\Base" FILES ${HEADER_FILE})
endforeach(HEADER_FILE)

