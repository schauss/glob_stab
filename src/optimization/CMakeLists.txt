include_directories(
    ../tools
    ../math
    ../parser
)

link_directories(
)

add_library( optimization SHARED
    function.cpp
    equation.cpp
    optimization.cpp
    optimization_result.cpp
    parameters.cpp
    value_set.cpp
)

target_link_libraries( optimization
    math
    parser
    tools
)

set(PROGRAMS )

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} optimization)
    set_target_properties(${program} PROPERTIES DEBUG_POSTFIX -debug)
    add_test(${program} ${program})
endforeach(program)

