POT Module
==========

**Namespace:** ``feff::pot``

**Header:** ``src/pot/pot.hpp``

The POT module calculates self-consistent muffin-tin potentials for each unique
atom type. This is Stage 2 of the pipeline.

Overview
--------

The module performs the following steps:

1. **Free-atom SCF** --- Solves the Dirac equation for each isolated atom type
   to obtain electron wave functions, densities, and energies.

2. **Muffin-tin construction** --- Overlaps the free-atom potentials using the
   Mattheiss prescription to construct the muffin-tin potential.

3. **Norman radius calculation** --- Determines the Norman (charge-neutral)
   radius for each atom type.

4. **SCF iteration** (optional) --- If the ``SCF`` card is present in the
   input, iterates the muffin-tin potentials to self-consistency using the
   full multiple scattering (FMS) Green's function.

Key Functions
-------------

``void pot(...)``
   Main computational kernel. Takes atomic coordinates, potential types,
   overlap parameters, and SCF flags. Produces muffin-tin radii, electron
   densities, and potentials. Writes ``pot.pad``.

Output
------

- ``pot.pad`` --- Binary file containing the converged potentials, densities,
  and muffin-tin geometry for use by XSPH.
- ``pot.json`` --- JSON metadata with SCF convergence information.
