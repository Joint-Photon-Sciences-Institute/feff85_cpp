add_library(fefffovrg
    nucdec.cpp
    radial_integrals_c.cpp
    intout.cpp
    solout.cpp
    solin.cpp
    wfirdc.cpp
    inmuac.cpp
    muatcc.cpp
    potex.cpp
    ortdac.cpp
    potdvp.cpp
    dfovrg.cpp
)

target_include_directories(fefffovrg
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(fefffovrg
    PUBLIC feffatom feffcommon feffmath feffpar
)
