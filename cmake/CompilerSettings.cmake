# Shared compiler flags for feff85exafs

if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
    # Static link C++ runtime so the executable is self-contained (no DLL dependencies)
    add_link_options(-static)
endif()

# Enable position-independent code for shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
