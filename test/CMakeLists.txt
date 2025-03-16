
if (ENABLE_TESTS)
    enable_testing()
    add_subdirectory(unit)
    add_subdirectory(integration)
    add_subdirectory(test_Utils)
endif()