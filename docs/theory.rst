Theory & Pipeline Architecture
==============================

What FEFF Calculates
--------------------

FEFF is an *ab initio* code for calculating X-ray Absorption Fine Structure
(XAFS), including both Extended X-ray Absorption Fine Structure (EXAFS) and
X-ray Absorption Near-Edge Structure (XANES). It computes the X-ray absorption
coefficient :math:`\mu(E)` and the EXAFS signal :math:`\chi(k)` using a
real-space multiple-scattering approach.

The EXAFS equation implemented by FEFF is:

.. math::

   \chi(k) = S_0^2 \sum_{\Gamma} \frac{N_\Gamma |f_\mathrm{eff}(k)|}{k R_\Gamma^2}
   \sin(2kR_\Gamma + \Phi_\Gamma(k)) \, e^{-2\sigma_\Gamma^2 k^2} \, e^{-2R_\Gamma/\lambda(k)}

where the sum runs over scattering paths :math:`\Gamma`, :math:`f_\mathrm{eff}`
is the effective curved-wave scattering amplitude, :math:`R_\Gamma` is the
effective half-path length, :math:`\sigma_\Gamma^2` is the Debye-Waller factor,
and :math:`\lambda(k)` is the photoelectron mean free path.

The 6-Stage Pipeline
--------------------

The FEFF calculation is organized as a sequential pipeline of six stages. Each
stage reads its input from JSON configuration files and binary pad files
produced by prior stages.

.. code-block:: text

   feff.inp
      |
      v
   +---------+     +-------+     +--------+     +--------+     +---------+     +-------+
   | 1.RDINP | --> | 2.POT | --> | 3.XSPH | --> | 4.PATH | --> | 5.GENFMT| --> | 6.FF2X|
   +---------+     +-------+     +--------+     +--------+     +---------+     +-------+
      |               |              |               |              |              |
   JSON configs    pot.pad       phase.pad       paths.dat      feff.pad      chi.dat
   atoms.json      pot.json      xsect.json      list.dat       feffNNNN.dat  xmu.dat
   global.json
   geom.json

Stage 1: RDINP --- Input Parsing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reads the ``feff.inp`` file and populates the ``FeffInput`` structure, which
contains all parameters organized by module. Writes JSON configuration files
(``atoms.json``, ``global.json``, ``geom.json``, ``pot.json``, etc.) that are
read by subsequent stages.

**Module:** :doc:`modules/rdinp`

Stage 2: POT --- Muffin-Tin Potentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculates the self-consistent muffin-tin potentials for each unique atom type.
This stage solves the atomic Dirac equation to obtain electron densities and
constructs overlapped muffin-tin potentials. If SCF is enabled, the potentials
are iterated to self-consistency.

Key computations:

- Atomic SCF for free atoms (Dirac equation)
- Muffin-tin radius determination
- Norman radius calculation
- Potential overlap (Mattheiss prescription)
- Optional self-consistent iteration

**Module:** :doc:`modules/pot`

Stage 3: XSPH --- Phase Shifts & Cross-Sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculates the scattering phase shifts :math:`\delta_l(E)` for each unique
potential type at each energy on the calculation grid. Also computes the
photoabsorption cross-section and the energy-dependent exchange-correlation
self-energy.

Key computations:

- Energy mesh construction (XANES + EXAFS regions)
- Phase shift calculation for each potential type and angular momentum
- Photoabsorption cross-section
- Full multiple scattering (FMS) for XANES region
- Writes ``phase.pad`` binary file

**Module:** :doc:`modules/xsph`

Stage 4: PATH --- Scattering Path Enumeration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Enumerates all significant multiple-scattering paths within a sphere of
radius ``RMAX`` around the absorbing atom. Uses a heap-based algorithm
with importance filtering to efficiently find paths that contribute
significantly to the EXAFS signal.

Key computations:

- Path enumeration using Rehr-Albers curved-wave theory
- Heap-based importance sampling
- Path degeneracy calculation
- Path filtering by curved-wave importance criterion

**Module:** :doc:`modules/path`

Stage 5: GENFMT --- Effective Scattering Amplitudes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculates the effective scattering amplitude (F-matrix) for each path found
in Stage 4. Uses the Rehr-Albers separable representation of the curved-wave
Green's function to compute the F-matrix efficiently.

Key computations:

- Rotation matrices for scattering geometry
- F-matrix calculation (Rehr-Albers approximation)
- Curved-wave corrections
- Writes ``feff.pad`` and individual ``feffNNNN.dat`` files

**Module:** :doc:`modules/genfmt`

Stage 6: FF2X --- Output Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Combines all path contributions to produce the final output spectra.
Sums the contributions from each scattering path, applying Debye-Waller
factors, mean free path corrections, and the amplitude reduction factor
:math:`S_0^2`.

Key computations:

- EXAFS :math:`\chi(k)` summation
- XANES :math:`\mu(E)` calculation
- Debye-Waller factor application (correlated Debye model)
- Self-energy corrections

**Outputs:** ``chi.dat``, ``xmu.dat``

**Module:** :doc:`modules/ff2x`

Data Flow
---------

The ``FeffInput`` structure (defined in ``include/feff/feff_input.hpp``) is the
central data container. It holds all parameters parsed from ``feff.inp`` and is
serialized to JSON files for inter-module communication. Each stage reads the
JSON files it needs and writes its results as binary pad files or additional
JSON files.

Physical Units
--------------

Internally, all calculations use **atomic units**:

- Distances in **Bohr radii** (1 Bohr = 0.52918 Angstrom)
- Energies in **Hartrees** (1 Hartree = 27.211 eV)
- Wave vectors in **inverse Bohr radii**

Input and output files use conventional units:

- Distances in **Angstroms**
- Energies in **eV**
- Wave vectors :math:`k` in **inverse Angstroms** (1/A)

The conversion constants are defined in ``include/feff/constants.hpp``
(see :doc:`constants`).

Radial Grid
-----------

FEFF uses a Loucks logarithmic radial grid for all radial integrations:

.. math::

   x(j) = -8.8 + (j - 1) \times 0.05, \qquad r(j) = e^{x(j)}

The grid has ``nrptx = 1251`` points, covering radii from
:math:`r \approx 1.5 \times 10^{-4}` Bohr to :math:`r \approx 86.9` Bohr.
This grid is used throughout the atomic SCF, potential construction, and phase
shift calculations.
