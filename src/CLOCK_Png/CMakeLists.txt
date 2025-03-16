

#Configure the install
 
find_package(LAPACK REQUIRED)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 

# include(CMakePrintHelpers)
# cmake_print_variables(CMAKE_INSTALL_PREFIX)

#Split png handling code into its own lib
add_library(Clock_Png STATIC
    Lib_Png.F90
    Lib_ColourScale.f90 
    Lib_Greyscale.c
)
set_target_properties(Clock_Png PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(Clock_Png
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(Clock_Png ${PNG_LIBRARY})

#Copy libraries to the lib directory of
#the install directory

install(TARGETS Clock_Png
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


 