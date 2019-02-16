include_directories(
    ../math
    ../parser
    ../optimization
    ../tools
)

link_directories(
)

set(PROGRAMS glob_stab)

find_package(Boost COMPONENTS system REQUIRED)

foreach(program ${PROGRAMS})
    add_executable(${program} ${program}.cpp)
	target_link_libraries(${program} math optimization parser tools ${Boost_SYSTEM_LIBRARY})
    set_target_properties(${program} PROPERTIES DEBUG_POSTFIX -debug)
endforeach(program)

