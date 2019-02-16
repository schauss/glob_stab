find_package( Boost COMPONENTS program_options thread)

include_directories(
)

add_library(tools SHARED
    branch_and_bound.hpp
    helpers.cpp
    config.cpp
    mat.hpp
    nd_array.hpp
    thread_pool.hpp
    vec.hpp
)

target_link_libraries(tools
    ${Boost_LIBRARIES}
)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)

set(PROGRAMS nd_array_test)

foreach(program ${PROGRAMS})
        add_executable(${program} ${program}.cpp)
        target_link_libraries(${program} tools ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
	set_target_properties(${program} PROPERTIES DEBUG_POSTFIX -debug)
        add_test(${program} ${program})
endforeach(program)
