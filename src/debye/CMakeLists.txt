add_library(feffdebye
    sigms.cpp
    sigcl.cpp
    sigm3.cpp
    sigte3.cpp
    sigrem.cpp
)

target_include_directories(feffdebye
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffdebye
    PUBLIC feffcommon feffmath
)
