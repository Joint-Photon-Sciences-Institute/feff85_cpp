add_library(feffmath
    distance.cpp
    bessel.cpp
    legendre.cpp
    wigner.cpp
    sommerfeld.cpp
    polynomial_roots.cpp
    determinant.cpp
    convolution.cpp
    interpolation.cpp
    integration.cpp
    phase_amplitude.cpp
    bcoef.cpp
)

target_include_directories(feffmath PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(feffmath PUBLIC
    Eigen3::Eigen
)
