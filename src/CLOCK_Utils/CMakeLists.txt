

#Configure the install
 
find_package(LAPACK REQUIRED)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 

# include(CMakePrintHelpers)
# cmake_print_variables(CMAKE_INSTALL_PREFIX)

#Split utility code into its own lib
add_library(Clock_Utils STATIC
    NBAX_StringTokenizers.f90
    Lib_CommandLineArguments.f90
    Lib_SimpleProgressBar.f90
    Lib_Callipers.f90
    Lib_Quicksort.f90
    Lib_RandomSeed.f90
    Lib_Filenames.f90
    Lib_ColouredTerminal.f90)
set_target_properties(Clock_Utils PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(Clock_Utils
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

#Copy libraries to the lib directory of
#the install directory

install(TARGETS Clock_Utils
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


 