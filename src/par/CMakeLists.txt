add_library(feffpar
    parallel.cpp
)

target_include_directories(feffpar
    PUBLIC ${CMAKE_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

option(FEFF_ENABLE_MPI "Enable MPI support" OFF)
if(FEFF_ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(feffpar PUBLIC MPI::MPI_CXX)
    target_compile_definitions(feffpar PUBLIC FEFF_USE_MPI)
endif()
