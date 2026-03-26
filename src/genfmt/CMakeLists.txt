add_library(feffgenfmt
    genfmt.cpp
    onepath.cpp
    genfmt_prep.cpp
    rdpath.cpp
    setlam.cpp
    rot3i.cpp
    mmtr.cpp
    mmtrxi.cpp
    sclmz.cpp
    fmtrxi.cpp
    snlm.cpp
    import_calc.cpp
    regenf.cpp
    xstar.cpp
    pathgeom.cpp
    read_onepath.cpp
    trig_utils.cpp
)

target_include_directories(feffgenfmt
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffgenfmt
    PUBLIC feffcommon feffmath feffpar feffjson
)
