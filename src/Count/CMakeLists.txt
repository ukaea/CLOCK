add_executable(count  analyseSpots.f90 countSpots.f90)
target_link_libraries(count Clock_Utils Clock_Png Clock_Img Clock fftw3 ${PNG_LIBRARY} ${LAPACK_LIBRARIES})

#Install the exec to bin folder of install path
install(TARGETS count
    RUNTIME
    DESTINATION bin)
    
