add_library(feffff2x
    reff2x.cpp
    rdfbin.cpp
    rdxbin.cpp
    exconv.cpp
    xscorr.cpp
    fprime.cpp
    feffdt.cpp
    dwadd.cpp
    ff2chi.cpp
    ff2xmu.cpp
    ff2afs.cpp
    ff2x.cpp
)

target_include_directories(feffff2x
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffff2x
    PUBLIC feffcommon feffmath feffdebye feffjson nlohmann_json::nlohmann_json
)
