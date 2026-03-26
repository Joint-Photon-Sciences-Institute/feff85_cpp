add_executable(math_tests
    math_test.cpp
)

target_link_libraries(math_tests
    feffmath
    GTest::gtest_main
)

add_executable(common_tests
    common_test.cpp
)

target_link_libraries(common_tests
    feffcommon
    feffpar
    feffmath
    GTest::gtest_main
)

add_executable(par_tests
    par_test.cpp
)

target_link_libraries(par_tests
    feffpar
    GTest::gtest_main
)

add_executable(atom_tests
    atom_test.cpp
)

target_link_libraries(atom_tests
    feffatom
    feffcommon
    feffmath
    feffpar
    GTest::gtest_main
)

add_executable(fovrg_tests
    fovrg_test.cpp
)

target_link_libraries(fovrg_tests
    fefffovrg
    feffatom
    feffcommon
    feffmath
    feffpar
    GTest::gtest_main
)

add_executable(exch_tests
    exch_test.cpp
)

target_link_libraries(exch_tests
    feffexch
    feffcommon
    feffmath
    feffpar
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(math_tests)
gtest_discover_tests(common_tests)
gtest_discover_tests(par_tests)
gtest_discover_tests(atom_tests)
gtest_discover_tests(fovrg_tests)
gtest_discover_tests(exch_tests)
