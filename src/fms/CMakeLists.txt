add_library(fefffms
    xprep.cpp
    yprep.cpp
    fms_core.cpp
    gglu.cpp
    ggbi.cpp
    ggrm.cpp
    gggm.cpp
    ggtf.cpp
    fmsie.cpp
    fmstot.cpp
    reafms.cpp
)

target_include_directories(fefffms
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(fefffms
    PUBLIC feffcommon feffmath feffpar Eigen3::Eigen nlohmann_json::nlohmann_json
)
