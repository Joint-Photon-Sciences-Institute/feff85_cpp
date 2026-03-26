FOVRG Module
============

**Namespace:** ``feff::fovrg``

**Header:** ``src/fovrg/dfovrg.hpp``

The FOVRG (Fine-structure OVerlay ReGion) module solves the Dirac equation in
the overlay region of the muffin-tin potential to obtain the scattering phase
shifts at each energy.

Overview
--------

This module solves the radial Dirac equation for a given muffin-tin potential
and energy, integrating from the origin outward and from the boundary inward,
then matching at an intermediate point. It handles both real (bound-state) and
complex (scattering) energies.

Key Functions
-------------

``void dfovrg(...)``
   Main Dirac equation solver. Integrates the coupled first-order differential
   equations for the large and small components of the Dirac spinor.

``void solout(...)``
   Outward integration of the Dirac equation from the origin.

``void intout(...)``
   Inward integration of the Dirac equation from the boundary.

Data Structures
---------------

``FovrgState``
   Working state for the overlay region solver, containing:

   - ``DiracWorkspaceComplex`` --- complex-valued working arrays
   - ``MeshParamsComplex`` --- radial mesh with complex potential
