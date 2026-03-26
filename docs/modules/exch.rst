EXCH Module
===========

**Namespace:** ``feff::exch``

**Header:** ``src/exch/xcpot.hpp``

The EXCH module calculates exchange-correlation potentials and self-energies
for the muffin-tin potentials.

Overview
--------

FEFF supports several exchange-correlation models, selected by the ``EXCHANGE``
card in the input:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Index
     - Model
     - Description
   * - 0
     - Hedin-Lundqvist (HL)
     - Complex self-energy with plasmon-pole model (default)
   * - 1
     - Dirac-Hara (DH)
     - Real, energy-dependent exchange
   * - 2
     - Ground-state
     - Real, ground-state exchange-correlation
   * - 5
     - DH + HL correlation
     - Dirac-Hara exchange with HL correlation

Key Functions
-------------

``void xcpot(...)``
   Calculate the exchange-correlation potential for a given electron density
   and energy. Handles all exchange models.

``void mpse(...)``
   Calculate the many-pole self-energy for the Hedin-Lundqvist model.
