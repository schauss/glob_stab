include_directories(
    ../math
    ../optimization
    ../tools
)

link_directories(
)

add_library( parser SHARED
    parse_equation.cpp
    parse_identifier.cpp
    parse_parameter.cpp
)

target_link_libraries( parser
    math
    tools
)

set(PROGRAMS )

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
    target_link_libraries(${program} parser)
    set_target_properties(${program} PROPERTIES DEBUG_POSTFIX -debug)
    add_test(${program} ${program})
endforeach(program)

