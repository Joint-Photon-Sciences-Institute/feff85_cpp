Installation
============

Prerequisites
-------------

**To build from source:**

- CMake 3.20 or later
- C++17 compiler (GCC/MinGW recommended on Windows, GCC/Clang on Linux/macOS)
- Ninja or Make build system

**To run the benchmark tests:**

- Python 3 with ``numpy`` and ``matplotlib``

Windows (MSYS2/MinGW)
---------------------

1. Install `MSYS2 <https://www.msys2.org/>`_ and open the **MINGW64** shell.

2. Install the required packages:

   .. code-block:: bash

      pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja

3. Clone and build:

   .. code-block:: bash

      git clone https://github.com/Joint-Photon-Sciences-Institute/feff85_cpp.git
      cd feff85_cpp

      cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
      cmake --build build

4. The executable will be at ``build/apps/feff8l.exe``.

.. note::

   You must build from an MSYS2 MINGW64 shell (not CMD or PowerShell) so the
   compiler can find its runtime libraries during the build process.

Linux / macOS
-------------

.. code-block:: bash

   # Install dependencies (Ubuntu/Debian example)
   sudo apt install build-essential cmake ninja-build

   # Clone and build
   git clone https://github.com/Joint-Photon-Sciences-Institute/feff85_cpp.git
   cd feff85_cpp

   cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
   cmake --build build

The executable will be at ``build/apps/feff8l``.

Static Linking
--------------

The build statically links the C++ runtime on GCC/MinGW, so the resulting
executable is fully self-contained with no external DLL dependencies. Only
Windows system DLLs (``ntdll.dll``, ``KERNEL32.DLL``, ``msvcrt.dll``) are
required at runtime.

Pre-built Binaries
------------------

If a pre-built binary is provided, simply place ``feff8l`` (or ``feff8l.exe``)
somewhere on your PATH. No additional libraries are required.
