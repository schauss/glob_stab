#Packages
include_directories(
    ../tools
)

link_directories(
)

add_library(math SHARED
    bernstein.cpp
    bernstein_exponent.cpp
    bernstein_2d.cpp
    bernstein_2d_patch.cpp
    iv.cpp
    iv_sse2.cpp
    iv_sse2_rounding.hpp
    monomial.cpp
    polygon.cpp
    iv_polygon.cpp
    polynomial.cpp
    taylor_model.cpp
)

target_link_libraries(math
    CGAL
    tools
    ${Boost_SYSTEM_LIBRARY}
    gmp
    mpfr
)

find_package(Boost COMPONENTS unit_test_framework system REQUIRED)

set(PROGRAMS iv_test polygon_test polynomial_test taylor_model_test iv_polygon_test)

foreach(program ${PROGRAMS})
        add_executable(${program} ${program}.cpp)
		target_link_libraries(${program} math ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_SYSTEM_LIBRARY})
	set_target_properties(${program} PROPERTIES DEBUG_POSTFIX -debug)
        add_test(${program} ${program})
endforeach(program)

