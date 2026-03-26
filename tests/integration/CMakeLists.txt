# Integration tests: full pipeline comparison against Fortran baselines

# Path to the Fortran test/baseline directory
set(FORTRAN_TEST_DIR "${CMAKE_SOURCE_DIR}/../feff85exafs-fortran/tests"
    CACHE PATH "Path to Fortran baseline test directory")

# Working directories for integration test output
set(INTEGRATION_WORK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/work"
    CACHE PATH "Working directory for integration test output")

# ---- RDINP parsing tests ----
add_executable(rdinp_integ_tests
    rdinp_test.cpp
)

target_include_directories(rdinp_integ_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(rdinp_integ_tests
    feffrdinp
    feffjson
    feffcommon
    feffpar
    feffmath
    GTest::gtest_main
)

target_compile_definitions(rdinp_integ_tests PRIVATE
    FORTRAN_TEST_DIR="${FORTRAN_TEST_DIR}"
    INTEGRATION_WORK_DIR="${INTEGRATION_WORK_DIR}"
)

# ---- Full pipeline integration tests ----
add_executable(pipeline_integ_tests
    pipeline_test.cpp
)

target_include_directories(pipeline_integ_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(pipeline_integ_tests
    feffpipeline
    GTest::gtest_main
)

target_compile_definitions(pipeline_integ_tests PRIVATE
    FORTRAN_TEST_DIR="${FORTRAN_TEST_DIR}"
    INTEGRATION_WORK_DIR="${INTEGRATION_WORK_DIR}"
)

# ---- Register with CTest ----
include(GoogleTest)
gtest_discover_tests(rdinp_integ_tests
    PROPERTIES LABELS "integration"
)
gtest_discover_tests(pipeline_integ_tests
    PROPERTIES LABELS "integration"
    # Pipeline tests may take longer
    DISCOVERY_TIMEOUT 120
)
