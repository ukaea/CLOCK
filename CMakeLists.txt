#Change this to whatever cmake version you have
cmake_minimum_required(VERSION 3.14) 


#Change the name and the version as you like
project(CLOCK
    VERSION 4.3.0
    #Add C, CXX etc if interface moves to the same project
    LANGUAGES C Fortran)

    

# Grab Python, 3.0 or newer
#find_package(Python 3.0 REQUIRED
  #COMPONENTS Interpreter Development.Module NumPy)

  
  

# Grab the variables from a local Python installation
# F2PY headers
# execute_process(
#   COMMAND "${Python_EXECUTABLE}"
#   -c "import numpy.f2py; print(numpy.f2py.get_include())"
#   OUTPUT_VARIABLE F2PY_INCLUDE_DIR
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )
     

include(CMakePrintHelpers)
include(CTest)
cmake_print_variables(PROJECT_NAME)
cmake_print_variables(CMAKE_INSTALL_PREFIX)
cmake_print_variables(PROJECT_SOURCE_DIR)

# cmake_print_variables(Python_INCLUDE_DIRS)
# cmake_print_variables(F2PY_INCLUDE_DIR)
# cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

        

#Add the option to turn tests off i.e. cmake -DENABLE_TESTS=OFF ....
option(ENABLE_TESTS "Build unit tests" ON)

#   add a preprocessor macro so code knows when it has been compiled with -DCMAKE_BUILD_TYPE=Debug
if(CMAKE_BUILD_TYPE STREQUAL "DBG_verbose")
    message(STATUS "setting preprocessor macro DEBUG")
    add_definitions(-DDEBUG)
endif()

#   add a preprocessor macro so code knows to compile against openMPI
if ( CMAKE_SERIAL )
    message ("forcing serial compilation")
    set(NPROCS_TEST "1")
else()
    message ("parallel compilation")
    #   Grab MPI?
    find_package(MPI)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})      
    if (MPI_Fortran_HAVE_F08_MODULE)
        message(STATUS "MPI F08-compliant compiler found, will add -DMPI flag to compilation") 
        add_definitions(-DMPI)
        if(DEFINED ENV{NPROCS})
            set(NPROCS_TEST "$ENV{NPROCS}")
        endif()
        message(STATUS "setting ctest processor count to ${NPROCS_TEST}")
    else()
        message(WARNING "MPI F08-compliant compiler not found")    
    endif()
endif()


#Don't know if the flags work for ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    #Set some compiler flags
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-256")
    add_compile_options(
        #Added regardless of build options
        "-I/usr/include"
        #For CMAKE_BUILD_TYPE=Release
        "$<$<CONFIG:Release>:-O2;-march=native;>"
        #For general profiling (not gprof), prevents inlining, and adds debug info but still optimised
        "$<$<CONFIG:Perf>:-Og;-g;-pg;>"
        #For gprof, prevents inlining, and adds debug info but still optimised
        "$<$<CONFIG:Gprof>:-Og;-fno-inline;-fno-omit-frame-pointer;-g>"
        #debug
        "$<$<CONFIG:Debug>:-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow;-g;-Wall;>"
        #light debug (for debugging with out the DEBUG preprocessor macro)
        "$<$<CONFIG:DBG_verbose>:-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow;-g;-Wall;>"
        "$<$<CONFIG:Coverage>:-coverage>"
        )
    add_link_options( 

         "$<$<CONFIG:Debug>:-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow;-g;>"
         "$<$<CONFIG:Perf>:-pg;>"
         "$<$<CONFIG:Gprof>:-p;-pg>"
         "$<$<CONFIG:DBG_verbose>:-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow;-g;-Wall;>"
         "$<$<CONFIG:Coverage>:-coverage>"
         )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")    
    #Add flags for intel (am assuming the old and LLVM compilers have the same flags 
    #I think can use STREQUAL "Intel" and STREQUAL "IntelLLVM" to distinguish
    #if needed
endif()

#Normally you want to add -DCMAKE_INSTALL_PREFIX=/path/to/install/location
#but this will default to the build directory for convenience
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH "..." FORCE)
endif()

add_subdirectory(src)
#add_subdirectory(bin)

if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()


#This stuff is boiler plate to export the details of
#the project to the install directory
#once done another cmake project can
#easily link to the library
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION lib/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE
    ${PROJECT_NAME}::
    DESTINATION
    lib/cmake)
export(TARGETS Clock_Utils Clock_Png Clock_Img Clock FILE ${PROJECT_NAME}Config.cmake)
