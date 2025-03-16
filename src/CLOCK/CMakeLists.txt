#Configure the install
 
find_package(LAPACK REQUIRED)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

#Split Clock processing code into its own lib
add_library(Clock STATIC
    Lib_FFTW3f.F90
    Lib_Splines1.f90 
    Lib_DataDoubler.f90
    Lib_Maxima2d.F90
    Lib_Gaussian2d.F90
    Lib_MultipleGaussians2d.F90 
    Lib_FitGaussians1d.f90
    Lib_GoldenSection.F90 
    Lib_Kernels.f90
    Lib_Stats.f90
    Lib_Graphs.f90
)

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW fftw3 IMPORTED_TARGET)
if(FFTW_FOUND)
  target_link_libraries(Clock INTERFACE PkgConfig::FFTW
  Clock PRIVATE Clock_Img
  )
  message("      PkgConfig has found FFTW3")
else()  # some systems do not have fftw3 registered in pkg-config, in this case it falls back to good old -lfftw3 which can fail later
  if(NOT DEFINED FFTW_LIBRARIES)
    message(WARNING "FFTW not detected by pkg-config, using fftw3 or define -DFFTW_LIBRARIES (e.g. =-lfftw3) and -DFFTW_INCLUDE_DIRS (e.g. =/usr/include)")
    target_link_libraries(Clock INTERFACE fftw3
    Clock PRIVATE Clock_Img)
    target_include_directories(Clock INTERFACE ${FFTW_INCLUDE_DIRS})
  else() #<------------------------------------------------------------------I am not sure if this part is needed, remove?
    target_link_libraries(Clock INTERFACE ${FFTW_LIBRARIES}
    Clock PRIVATE Clock_Img) 
    target_include_directories(Clock INTERFACE ${FFTW_INCLUDE_DIRS})
  endif()
endif()

set_target_properties(Clock PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(Clock
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

#Copy libraries to the lib directory of
#the install directory

install(TARGETS Clock
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE
    DESTINATION lib
    INCLUDES 
    DESTINATION
    include)
#...and the fortran module to the include directory
install(DIRECTORY
    ${CMAKE_BINARY_DIR}/include
    DESTINATION
    ${CMAKE_INSTALL_PREFIX})






 