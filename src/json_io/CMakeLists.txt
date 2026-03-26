add_library(feffjson
    json_io.cpp
)

target_include_directories(feffjson
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffjson
    PUBLIC feffcommon feffmath feffpar nlohmann_json::nlohmann_json
)
