# C++ ATOM module sources — pure C++ implementation.

set(ATOM_SOURCES
    utility.cpp
    nucdev.cpp
    radial_integrals.cpp
    muatco.cpp
    intdir.cpp
    soldir.cpp
    wfirdf.cpp
    inmuat.cpp
    ortdat.cpp
    lagdat.cpp
    potrdf.cpp
    vlda.cpp
    etotal.cpp
    tabrat.cpp
    fpf0.cpp
    s02at.cpp
    scfdat.cpp
)

add_library(feffatom ${ATOM_SOURCES})

target_include_directories(feffatom
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffatom
    PUBLIC feffcommon feffmath feffpar
)
