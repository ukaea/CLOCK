add_executable(imageCorrelation imageCorrelation.f90)
target_link_libraries(imageCorrelation Clock_Utils Clock_Png Clock_Img ${PNG_LIBRARY})

#Install the exec to bin folder of install path
install(TARGETS imageCorrelation
    RUNTIME
    DESTINATION bin)
    
