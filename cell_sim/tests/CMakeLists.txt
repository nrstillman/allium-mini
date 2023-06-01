add_executable(
    simulation_tests
    simulation_tests.cpp
)

target_link_libraries(
    simulation_tests
    LINK_PUBLIC
    capmd
    ${VTK_LIBRARIES}
)
