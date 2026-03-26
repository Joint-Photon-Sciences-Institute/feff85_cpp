add_library(feffpipeline
    feff_pipeline.cpp
)

target_include_directories(feffpipeline
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(feffpipeline
    PUBLIC
    feffff2x
    feffgenfmt
    feffpath
    feffxsph
    feffpot
    fefffms
    feffpot
    feffexch
    fefffovrg
    feffatom
    feffrdinp
    feffjson
    feffdebye
    feffcommon
    feffpar
    feffmath
)
