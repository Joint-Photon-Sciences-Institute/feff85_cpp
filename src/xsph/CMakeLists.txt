add_library(feffxsph
    xsph.cpp
    phase.cpp
    xsect.cpp
    phmesh.cpp
    phmesh2.cpp
    wrxsph.cpp
    rexsph.cpp
    axafs.cpp
    acoef.cpp
    radint.cpp
    szlz.cpp
    rholat.cpp
    rholsz.cpp
    fmssz.cpp
    wphase.cpp
    getedg.cpp
    rdgrid.cpp
)

target_include_directories(feffxsph
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffxsph
    PUBLIC fefffms fefffovrg feffexch feffjson feffcommon feffmath feffpar nlohmann_json::nlohmann_json
)
