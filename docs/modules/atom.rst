ATOM Module
===========

**Namespace:** ``feff::atom``

**Header:** ``src/atom/atom_types.hpp``

The ATOM module performs self-consistent field (SCF) calculations for free atoms
by solving the relativistic Dirac equation. It provides the starting electron
densities and potentials used by the POT module.

Overview
--------

For each unique atom type in the calculation, the ATOM module:

1. Sets up the electron configuration (shells, occupations, quantum numbers)
2. Solves the Dirac equation on a logarithmic radial grid
3. Iterates to self-consistency
4. Produces radial wave functions, electron densities, and total energies

Key Data Structures
-------------------

``AtomState``
   Aggregate state for an atomic SCF calculation, containing:

   - ``OrbitalArraysReal`` --- radial wave functions (large/small components)
   - ``OrbitalConfig`` --- electron configuration and energies
   - ``ScfParams`` --- SCF convergence parameters
   - ``DiracWorkspaceReal`` --- working arrays for Dirac solver
   - ``MeshParamsReal`` --- radial mesh parameters
   - ``LagrangeParams`` --- Lagrange multiplier data

``OrbitalArraysReal``
   Wave functions on a 251-point grid: ``cg[251][30]`` (large component),
   ``cp[251][30]`` (small component).

``OrbitalConfig``
   Quantum numbers (``nqn``, ``nk``, ``nel``), occupation numbers, and
   orbital energies for up to 30 orbitals.

Key Functions
-------------

``void scfdat(...)``
   Main atomic SCF driver. Iterates the Dirac equation to self-consistency.

``void etotal(...)``
   Compute the total atomic energy.

``void getorb(...)``
   Look up the electron configuration for a given element.

``void muatco(...)``
   Calculate the atomic absorption coefficient.

``void inmuat(...)``
   Initialize atomic data structures.
