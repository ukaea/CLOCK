
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
 


#Here is a trivial example test that just checks libpng still works
add_executable(testpng testPng.f90)
target_link_libraries(testpng test_Utils Clock_Utils Clock_Png ${PNG_LIBRARY} )
add_test(NAME pngtest COMMAND testpng -f ${PROJECT_SOURCE_DIR}/data/test_cameraman.png )
set_tests_properties(pngtest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
 
add_executable(testrandc testRandC.f90)
target_link_libraries(testrandc test_Utils Clock_Utils Clock_Img )
add_test(NAME randctest COMMAND testrandc)
set_tests_properties(randctest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
     
    
     
find_program(BASH_PROGRAM bash)
if (BASH_PROGRAM)
  add_test(counttest ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/testCount.sh ${PROJECT_SOURCE_DIR}/data/)
endif()
set_tests_properties(counttest
    PROPERTIES
    LABELS "EXAMPLE;LONG"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 


find_program(BASH_PROGRAM bash)
if (BASH_PROGRAM)
  add_test(bfdftest ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/testBfdf.sh ${PROJECT_SOURCE_DIR}/data/)
endif()
set_tests_properties(bfdftest
    PROPERTIES
    LABELS "EXAMPLE;LONG"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")    
    

add_executable(testFilenames testLib_Filenames.f90)
target_link_libraries(testFilenames test_Utils Clock_Utils)
add_test(NAME Filenamestest COMMAND testFilenames )
set_tests_properties(Filenamestest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")  


add_executable(testQuicksort testLib_Quicksort.f90)
target_link_libraries(testQuicksort test_Utils Clock_Utils )
add_test(NAME Quicksorttest COMMAND testQuicksort)
set_tests_properties(Quicksorttest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    ##ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")


add_executable(testTokenizer testNBAX_StringTokenizers.f90)
target_link_libraries(testTokenizer test_Utils Clock_Utils)
add_test(NAME Tokenizertest COMMAND testTokenizer)
set_tests_properties(Tokenizertest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testGoldenSelection testLib_GoldenSection.f90)
target_link_libraries(testGoldenSelection test_Utils Clock_Utils Clock)
add_test(NAME GoldenSelectiontest COMMAND testGoldenSelection)
set_tests_properties(GoldenSelectiontest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")


add_executable(testRadon testRadon.f90)
target_link_libraries(testRadon test_Utils Clock_Utils Clock_Img Clock_Png)
add_test(NAME Radontest COMMAND testRadon ${PROJECT_SOURCE_DIR}/data/test_line.png)
set_tests_properties(Radontest
    PROPERTIES
    LABELS "EXAMPLE;LONG"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")
    

add_executable(testFilter testFilter.f90)
target_link_libraries(testFilter test_Utils Clock_Utils Clock_Png Clock_Img Clock fftw3 ${PNG_LIBRARY} ${LAPACK_LIBRARIES})
add_test(NAME FilterTest COMMAND testFilter -f ${PROJECT_SOURCE_DIR}/data/filtertest.png -o ${PROJECT_SOURCE_DIR}/data/filtertest.sandpTEST.png -sandp -lambda 2.0 -noflat -hpf 2)
set_tests_properties(FilterTest
    PROPERTIES
    LABELS "EXAMPLE;LONG"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 

add_executable(testFlatten testFlatten.f90)
target_link_libraries(testFlatten test_Utils Clock_Utils Clock_Png Clock_Img Clock fftw3 ${PNG_LIBRARY} ${LAPACK_LIBRARIES})
add_test(NAME FlattenTest COMMAND testFlatten -f ${PROJECT_SOURCE_DIR}/data/flattest.png -o ${PROJECT_SOURCE_DIR}/data/flattest.flatTEST.png -nozero -n 4)
set_tests_properties(FlattenTest
    PROPERTIES
    LABELS "EXAMPLE;LONG"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(test_saltandpepperfilter test_saltandpepperfilter.f90)
target_link_libraries(test_saltandpepperfilter Clock_Utils Clock_Img )
add_test( NAME saltandpepperfilter_test1 COMMAND test_saltandpepperfilter -t sandp1 -f ${CMAKE_INSTALL_PREFIX}/src/SaltAndPepperFilter/saltAndPepperFilter)
set_tests_properties(saltandpepperfilter_test1
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
add_test( NAME saltandpepperfilter_test3 COMMAND test_saltandpepperfilter -t sandp3 -f ${CMAKE_INSTALL_PREFIX}/src/SaltAndPepperFilter/saltAndPepperFilter)
set_tests_properties(saltandpepperfilter_test3
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 

add_executable(test_mlefilter test_MLEFilter.f90)
target_link_libraries(test_mlefilter test_Utils Clock_Utils Clock_Png Clock_Img  ${PNG_LIBRARY})
add_test( NAME mlefiltertest COMMAND test_mlefilter)
set_tests_properties(mlefiltertest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS") 
 