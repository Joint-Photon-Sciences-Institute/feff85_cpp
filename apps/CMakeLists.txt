# Standalone executables

add_executable(feff8l feff8l.cpp)
target_link_libraries(feff8l PRIVATE feffpipeline)

# Increase stack size for large arrays in POT/ATOM calculations
if(MSVC)
    target_link_options(feff8l PRIVATE /STACK:16777216)
else()
    target_link_options(feff8l PRIVATE -Wl,--stack,67108864)
endif()

# Individual module drivers (optional)
# add_executable(feff_pot feff_pot.cpp)
# target_link_libraries(feff_pot PRIVATE feffpipeline)
