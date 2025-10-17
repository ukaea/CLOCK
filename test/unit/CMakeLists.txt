
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

     
find_program(BASH_PROGRAM bash)
if (BASH_PROGRAM)
  add_test(runbin ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/canRunBin.sh )
endif()
set_tests_properties(runbin
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")     

 
add_executable(testCLA testLib_CommandLineArguments.f90)
target_link_libraries(testCLA test_Utils Clock_Utils )
add_test(NAME CLAtest COMMAND testCLA -a 12345 -b 1,2,3,4,5 -c 1.2345E+06 -d 1.0e-5,-0.0,0.0,12345.6789,1.0e+5 -e -f "sandwiches, herring")
set_tests_properties(CLAtest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")    


add_executable(testcallipers testLib_Callipers.f90)
target_link_libraries(testcallipers test_Utils Clock_Utils )
add_test(NAME calliperstest COMMAND testcallipers)
set_tests_properties(calliperstest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testRandom testLib_RandomSeed.f90)
target_link_libraries(testRandom test_Utils Clock_Utils )
add_test(NAME Randomtest COMMAND testRandom)
set_tests_properties(Randomtest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")


add_executable(testColourTerminal testLib_ColouredTerminal.f90)
target_link_libraries(testColourTerminal test_Utils Clock_Utils)
add_test(NAME ColourTerminalTest COMMAND testColourTerminal)
set_tests_properties(ColourTerminalTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testColourScale testLib_ColourScale.f90)
target_link_libraries(testColourScale test_Utils Clock_Png Clock_Utils)
add_test(NAME ColourScaleTest COMMAND testColourScale)
set_tests_properties(ColourScaleTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testStats testLib_Stats.f90)
target_link_libraries(testStats Clock test_Utils Clock_Utils)
add_test(NAME StatsTest COMMAND testStats)
set_tests_properties(StatsTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")


add_executable(testGraphs testLib_Graphs.f90)
target_link_libraries(testGraphs test_Utils Clock Clock_Utils)
add_test(NAME GraphsTest COMMAND testGraphs)
set_tests_properties(GraphsTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testImageTools testLib_ImageTools.f90)
target_link_libraries(testImageTools test_Utils Clock Clock_Img)
add_test(NAME ImageToolsTest COMMAND testImageTools)
set_tests_properties(ImageToolsTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testLib_SafeExp testLib_SafeExp.f90)
target_link_libraries(testLib_SafeExp Clock_Utils )
add_test(NAME testLib_SafeExp_test COMMAND testLib_SafeExp)
set_tests_properties(testLib_SafeExp_test
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testDrawEllipse testLib_DrawEllipse.f90)
target_link_libraries(testDrawEllipse test_Utils Clock Clock_Img)
add_test(NAME DrawEllipseTest COMMAND testDrawEllipse)
set_tests_properties(DrawEllipseTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

add_executable(testHistograms testLib_Histograms.f90)
target_link_libraries(testHistograms test_Utils Clock Clock_Img)
add_test(NAME HistogramsTest COMMAND testHistograms)
set_tests_properties(HistogramsTest
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS")

#CLEAN UP NEEDS TO BE LAST
add_test(NAME clean_up_unit_tests #this calls cleanuptests.cmake that deletes all temporary files associated with testing 
    COMMAND ${CMAKE_COMMAND} #this might need to be last
      -P ${CMAKE_CURRENT_SOURCE_DIR}/clean_up_unit_tests.cmake
  )