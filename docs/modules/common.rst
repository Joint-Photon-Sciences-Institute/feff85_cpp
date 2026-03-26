Common Module
=============

**Namespace:** ``feff::common``

**Headers:** ``src/common/``

The Common module provides shared utilities used throughout the codebase:
logging, file I/O, string manipulation, radial grid functions, and physics
helper routines.

Logging
-------

**Header:** ``src/common/logging.hpp``

``Logger& logger()``
   Returns the global logger instance. Writes to both ``stdout`` and
   ``feff.log``.

``void wlog(const std::string& msg)``
   Write a message to the log. Respects parallel execution settings
   (suppresses output on non-root processes).

Radial Grid
-----------

**Header:** ``src/common/radial_grid.hpp``

FEFF uses a Loucks logarithmic grid with ``nrptx = 1251`` points:

``double xx(int j)``
   Grid x-value: :math:`x(j) = -8.8 + (j-1) \times 0.05`

``double rr(int j)``
   Radial coordinate: :math:`r(j) = e^{x(j)}`

``int ii(double r)``
   Grid index from radius: :math:`j = 1 + (\\log r + 8.8) / 0.05`

File I/O
--------

**Header:** ``src/common/file_io.hpp``

Utilities for reading and writing FEFF's binary pad files and text data files.

**Header:** ``src/common/pad_io.hpp``

Functions for reading and writing the binary ``pot.pad``, ``phase.pad``, and
``feff.pad`` files.

String Utilities
----------------

**Header:** ``src/common/string_utils.hpp``

Helper functions for parsing ``feff.inp`` cards: trimming, case conversion,
tokenizing, and number parsing.

Physics Utilities
-----------------

**Header:** ``src/common/physics_utils.hpp``

Shared physics helper functions: Fermi function, electron density calculations,
and unit conversions.

Grid Interpolation
------------------

**Header:** ``src/common/grid_interpolation.hpp``

Interpolation routines for transferring data between different radial grids.

Orbital Data
------------

**Header:** ``src/common/orbital_data.hpp``

``void getorb(int iz, int ihole, double xion, int iunf, ...)``
   Look up the ground-state electron configuration for element ``iz`` with
   core hole ``ihole`` and ionization ``xion``. Returns quantum numbers,
   occupation numbers, and orbital data for all shells.
