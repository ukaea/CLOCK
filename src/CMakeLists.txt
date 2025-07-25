#Configure the install
 
#find_package(LAPACK REQUIRED)

#find_package(PNG REQUIRED)
#include_directories(${PNG_INCLUDE_DIR})
 

# include(CMakePrintHelpers)
# cmake_print_variables(CMAKE_INSTALL_PREFIX)

#Subdirectories for the executables
add_subdirectory(BFDF)
add_subdirectory(Count)
add_subdirectory(Flat)
add_subdirectory(Linear)
add_subdirectory(SaltAndPepperFilter)
add_subdirectory(MLE_Filter)

#Subdirectories for the static libraries
add_subdirectory(CLOCK_Utils)
add_subdirectory(CLOCK_Img)
add_subdirectory(CLOCK_Png)
add_subdirectory(CLOCK)
