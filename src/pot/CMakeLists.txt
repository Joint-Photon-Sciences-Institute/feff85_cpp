add_library(feffpot
    inipot.cpp
    moveh.cpp
    fermi.cpp
    frnrm.cpp
    sumax.cpp
    ovrlp.cpp
    sidx.cpp
    istprm.cpp
    afolp.cpp
    corval.cpp
    grids.cpp
    broydn.cpp
    coulom.cpp
    ff2g.cpp
    rholie.cpp
    wpot.cpp
    reapot.cpp
    scmt.cpp
    pot.cpp
    wrpot.cpp
    ovp2mt.cpp
    movrlp.cpp
    scmtmp.cpp
    lapack_shim.cpp
)

target_include_directories(feffpot
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffpot
    PUBLIC feffatom fefffovrg feffexch feffjson feffcommon feffmath feffpar
)
