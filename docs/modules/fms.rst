FMS Module
==========

**Namespace:** ``feff::fms``

**Header:** ``src/fms/fms_core.hpp``

The FMS (Full Multiple Scattering) module computes the Green's function for
a finite atomic cluster using matrix inversion. It is used in both the SCF
cycle (Stage 2) and the XANES calculation (Stage 3).

Overview
--------

FMS solves the multiple scattering equation:

.. math::

   G = G_0 + G_0 \, T \, G

by direct matrix inversion of :math:`(1 - G_0 T)` in the angular momentum
basis :math:`|i, l, m\rangle`, where :math:`i` labels atoms in the cluster
and :math:`l, m` are angular momentum quantum numbers.

The cluster is defined by a sphere of radius ``RFMS`` around the absorbing
atom. The angular momentum expansion is truncated at ``LFMS``.

Key Functions
-------------

``void fms(...)``
   Compute the FMS Green's function at a single energy. Returns the
   Green's function matrix in the angular momentum basis.

``void xprep(...)``
   Prepare the cluster geometry for FMS: select atoms within the FMS
   radius, compute interatomic distances and angles.

``void fmsie(...)``
   FMS matrix inversion engine. Constructs and inverts the
   :math:`(1 - G_0 T)` matrix.

``void getkts(...)``
   Build the basis states :math:`|i, l, m, s\rangle` for the cluster.

Data Structures
---------------

``FMSData``
   Aggregate shared state for FMS calculations:

   - ``ClusterData`` --- atom coordinates and potential types
   - ``RotationData`` --- rotation matrices for coordinate transformations
   - ``BasisStates`` --- quantum state labels
   - ``ClebschGordon`` --- Clebsch-Gordan coupling coefficients
   - ``DebyeWaller`` --- thermal displacement parameters

Key Parameters
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Parameter
     - Description
   * - ``rfms1``
     - FMS radius for SCF (Stage 2, typically 3-5 Angstroms)
   * - ``rfms2``
     - FMS radius for XANES (Stage 3, typically 5-8 Angstroms)
   * - ``lfms1``
     - Max angular momentum for SCF FMS
   * - ``lfms2``
     - Max angular momentum for XANES FMS (max 4)
