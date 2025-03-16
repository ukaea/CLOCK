add_executable(linear LinearFeatureDetect.f90 linearFeatures.f90)
target_link_libraries(linear Clock_Utils Clock_Png  Clock_Img Clock fftw3 ${PNG_LIBRARY} ${LAPACK_LIBRARIES})

#Install the exec to bin folder of install path
install(TARGETS linear
    RUNTIME
    DESTINATION bin)
    
