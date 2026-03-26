PATH Module
===========

**Namespace:** ``feff::path``

**Header:** ``src/path/paths.hpp``

The PATH module enumerates all significant multiple-scattering paths within a
cluster. This is Stage 4 of the pipeline.

Overview
--------

The path finder uses a heap-based algorithm to efficiently enumerate scattering
paths up to a maximum half-path length ``RMAX`` and a maximum number of legs
``NLEG``. Each path is characterized by:

- The sequence of scattering atoms visited
- The total path length :math:`2R`
- The degeneracy (number of equivalent paths)
- An importance criterion based on curved-wave scattering amplitudes

Paths are filtered by importance to retain only those that contribute
significantly to the EXAFS signal.

Key Functions
-------------

``void paths(...)``
   Main path enumeration routine. Uses a heap data structure to efficiently
   explore the space of possible scattering paths.

Key Parameters
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Parameter
     - Description
   * - ``rmax``
     - Maximum half-path length (Angstroms)
   * - ``nleg``
     - Maximum number of legs (scattering events + 1)
   * - ``pcritk``
     - Output filter: keep paths with importance > pcritk% of largest
   * - ``pcrith``
     - Heap filter: explore paths with importance > pcrith% of largest

Output
------

- ``paths.dat`` --- Detailed list of all significant paths with geometry.
- ``list.dat`` --- Compact path list with indices and degeneracies.
