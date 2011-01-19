if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    
  else(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
  
	find_library(MKL_intel_c
	  NAMES mkl_intel_c
	  PATHS "C:\Programme\Intel\MKL\9.0\ia32\lib" 
          "C:\Programme\Intel\MKL\9.0\ia32\bin"
          $ENV{LIB}
	)
  
    find_library(MKL_core
	  NAMES mkl_core
	  PATHS "C:\Programme\Intel\MKL\9.0\ia32\lib" 
          "C:\Programme\Intel\MKL\9.0\ia32\bin"
          $ENV{LIB}
	)
	
   find_library(MKL_intel_thread
	  NAMES mkl_intel_thread
	  PATHS "C:\Programme\Intel\MKL\9.0\ia32\lib" 
          "C:\Programme\Intel\MKL\9.0\ia32\bin"
          $ENV{LIB}
	)
	
    find_library(MKL_guide
	  NAMES libguide
	  PATHS "C:\Programme\Intel\MKL\9.0\ia32\lib" 
          "C:\Programme\Intel\MKL\9.0\ia32\bin"
          $ENV{LIB}
	)

	if(MKL_intel_c AND MKL_core AND MKL_intel_thread AND MKL_guide)
#	if(MKL_intel_c AND MKL_core AND MKL_guide)
	  message(STATUS "MKL Found")
	  set(MKL_FOUND TRUE)
	  set(MKL_LIBRARIES ${MKL_intel_c} ${MKL_core} ${MKL_intel_thread} ${MKL_guide})
#	  set(MKL_LIBRARIES ${MKL_intel_c} ${MKL_core} ${MKL_guide})
	  message(STATUS ${MKL_LIBRARIES})
      mark_as_advanced(MKL_LIBRARIES)
      mark_as_advanced(MKL_FOUND)
	endif(MKL_intel_c AND MKL_core AND MKL_intel_thread AND MKL_guide)
#	endif(MKL_intel_c AND MKL_core AND MKL_guide)
	  
    #include(FindPackageHandleStandardArgs)
    #find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES)

  endif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
  
else(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

endif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
