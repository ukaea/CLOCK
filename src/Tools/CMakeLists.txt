add_executable(scaleImage scaleImage.f90)
target_link_libraries(scaleImage Clock_Utils Clock_Png Clock ${PNG_LIBRARY})

add_executable(imageUtilities imageUtilities.f90)
target_link_libraries(imageUtilities Clock_Utils Clock_Png Clock_Img Clock ${PNG_LIBRARY})

#Install the exec to bin folder of install path
install(TARGETS scaleImage
    RUNTIME
    DESTINATION bin)
    
