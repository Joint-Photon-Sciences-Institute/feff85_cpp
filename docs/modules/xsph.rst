XSPH Module
===========

**Namespace:** ``feff::xsph``

**Header:** ``src/xsph/xsph.hpp``

The XSPH module calculates scattering phase shifts and the photoabsorption
cross-section. This is Stage 3 of the pipeline.

Overview
--------

1. **Energy mesh construction** --- Builds a combined energy grid with fine
   spacing in the XANES region and coarser spacing in the EXAFS region.

2. **Phase shift calculation** --- For each unique potential type and angular
   momentum :math:`l`, solves the radial Schrodinger/Dirac equation to obtain
   the scattering phase shifts :math:`\delta_l(E)`.

3. **Cross-section calculation** --- Computes the photoabsorption cross-section
   including matrix elements for the specified core hole.

4. **FMS for XANES** --- Applies full multiple scattering to obtain the XANES
   absorption using the cluster geometry.

Key Functions
-------------

``void xsph(...)``
   Main XSPH calculation. Takes muffin-tin potentials, atomic geometry, and
   exchange-correlation parameters. Produces phase shifts and cross-sections.

Supporting Functions
--------------------

``void xsect(...)``
   Calculate the photoabsorption cross-section for a single energy point.

``void phmesh2(...)``
   Construct the energy mesh for XANES and EXAFS regions.

``void fmssz(...)``
   Full multiple scattering in the spin-orbit basis.

``void szlz(...)``
   Spin-orbit coupling calculations.

``void rholat(...)``
   Radial integration for the lattice Green's function.

Output
------

- ``phase.pad`` --- Binary file containing phase shifts for all potential types.
- ``xsect.json`` --- Cross-section data and energy mesh metadata.
