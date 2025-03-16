add_executable(flat flattenBackground.f90)
target_link_libraries(flat Clock_Utils Clock_Png Clock_Img Clock fftw3 ${PNG_LIBRARY} ${LAPACK_LIBRARIES})

#Install the exec to bin folder of install path
install(TARGETS flat
    RUNTIME
    DESTINATION bin)
    
