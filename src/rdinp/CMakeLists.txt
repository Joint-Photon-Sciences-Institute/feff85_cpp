add_library(feffrdinp
    iniall.cpp
    setedg.cpp
    eu2au.cpp
    mkptz.cpp
    rdinp.cpp
    ffsort.cpp
)

target_include_directories(feffrdinp
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffrdinp
    PUBLIC feffcommon feffmath feffpar
)
