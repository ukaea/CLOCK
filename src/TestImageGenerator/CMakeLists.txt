
add_executable(testImageGenerator testImageGenerator.f90)
target_link_libraries(testImageGenerator Clock_Utils Clock_Png ${PNG_LIBRARY} Clock_Img Clock)

#Install the exec to bin folder of install path
install(TARGETS testImageGenerator 
    RUNTIME
    DESTINATION bin)
    
