add_library(feffexch
    vbh.cpp
    edp.cpp
    rhl.cpp
    imhl.cpp
    cubic_exch.cpp
    ffq.cpp
    quinn.cpp
    rhlbp.cpp
    fndsng.cpp
    csigma.cpp
    csigz.cpp
    xcpot.cpp
    mpse.cpp
)

target_include_directories(feffexch
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffexch
    PUBLIC feffcommon feffmath feffpar
)
