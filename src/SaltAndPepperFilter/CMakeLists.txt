
add_executable(saltAndPepperFilter saltAndPepperFilter.F90)
target_link_libraries(saltAndPepperFilter Clock_Utils Clock_Png ${PNG_LIBRARY})

#Install the exec to bin folder of install path
install(TARGETS saltAndPepperFilter 
    RUNTIME
    DESTINATION bin)
    
