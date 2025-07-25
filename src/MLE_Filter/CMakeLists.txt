
add_executable(MLEFilter MLEFilter.F90)
target_link_libraries(MLEFilter Clock_Utils Clock_Png Clock_Img ${PNG_LIBRARY})

#Install the exec to bin folder of install path
install(TARGETS MLEFilter 
    RUNTIME
    DESTINATION bin)
    
