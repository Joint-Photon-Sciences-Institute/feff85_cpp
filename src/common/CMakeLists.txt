add_library(feffcommon
    logging.cpp
    string_utils.cpp
    physics_utils.cpp
    orbital_data.cpp
    grid_interpolation.cpp
    pad_io.cpp
    file_io.cpp
    itoken.cpp
)

target_include_directories(feffcommon
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffcommon
    PUBLIC feffmath feffpar
)
