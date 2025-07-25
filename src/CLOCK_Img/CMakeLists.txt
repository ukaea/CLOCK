

#Configure the install
 
find_package(LAPACK REQUIRED)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 

# include(CMakePrintHelpers)
# cmake_print_variables(CMAKE_INSTALL_PREFIX)
   
#Split image processing code into its own lib
add_library(Clock_Img STATIC
    Lib_RidlerCalvard.F90
    Lib_MaxLikelihoodFilter.f90 
    Lib_LocalFilter2d.F90
    Lib_LowPassFilter3d.F90
    Lib_GaussianBlurs.F90
    Lib_Exteriors.F90 
    Lib_ImageCorrelationFunction.f90
    Lib_ConjugateGradient.F90
    Lib_LaplacianOfGaussians.F90    
    Lib_Radon.f90
    Lib_ImageTools.f90
    Lib_ExtrapolateImage.F90
    Lib_LaplacianKernel.f90

)
set_target_properties(Clock_Img PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(Clock_Img
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(Clock_Img Clock_Png Clock_Utils ${PNG_LIBRARY} ${LAPACK_LIBRARIES} )#need to sperate out Png stuff in Img back to Png

#Copy libraries to the lib directory of
#the install directory

install(TARGETS Clock_Img
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


 