GENFMT Module
=============

**Namespace:** ``feff::genfmt``

**Header:** ``src/genfmt/genfmt.hpp``

The GENFMT module calculates the effective scattering amplitudes (F-matrix)
for each path. This is Stage 5 of the pipeline.

Overview
--------

GENFMT implements the Rehr-Albers (RA) separable representation of curved-wave
multiple scattering. For each path, it computes the effective scattering
amplitude :math:`f_\mathrm{eff}(k)` and the total scattering phase
:math:`\Phi(k)` as a function of photoelectron wave vector :math:`k`.

The RA approximation factorizes the curved-wave propagator into products of
matrices at each scattering vertex, making the computation scale linearly
with the number of legs rather than exponentially.

Key Functions
-------------

``void genfmt(...)``
   Main F-matrix calculation for all paths. Reads ``phase.pad`` and
   ``paths.dat``, writes ``feff.pad`` and individual ``feffNNNN.dat`` files.

``void onepath(...)``
   Calculate the F-matrix for a single scattering path.

``void fmtrxi(...)``
   Compute the F-matrix elements using the RA separable approximation.

``void mmtr(...)``
   Compute the M-matrix (transfer matrix) for a single scattering vertex.

Data Structures
---------------

``PhaseData``
   Phase shifts and energy mesh read from ``phase.pad``.

``FmatrixData``
   F-matrix elements for a single path.

``RotationMatrixData``
   Rotation matrices (Wigner d-functions) for coordinate transformations
   between successive scattering planes.

``BmatrixData``
   B-matrix workspace for transfer matrix calculations.

Output
------

- ``feff.pad`` --- Binary file with combined path amplitudes and phases.
- ``feffNNNN.dat`` --- Individual path files with :math:`f_\mathrm{eff}(k)`,
  phase, and other path-specific data.
