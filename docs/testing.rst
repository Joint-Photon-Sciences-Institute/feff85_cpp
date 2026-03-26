Testing & Benchmarks
====================

The repository includes integration tests that compare the C++ output against
the original Fortran FEFF85 baseline for 8 materials.

Test Materials
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Material
     - Edge
     - Description
   * - Copper
     - Cu K
     - FCC metal
   * - NiO
     - Ni K
     - Nickel oxide
   * - UO2
     - U L3
     - Uranium dioxide
   * - LCO-para
     - Cu K
     - La\ :sub:`2`\ CuO\ :sub:`4`, parallel polarization
   * - LCO-perp
     - Cu K
     - La\ :sub:`2`\ CuO\ :sub:`4`, perpendicular polarization
   * - Zircon
     - Zr K
     - ZrSiO\ :sub:`4`
   * - bromoadamantane
     - Br K
     - C\ :sub:`10`\ H\ :sub:`15`\ Br (molecular)
   * - ferrocene
     - Fe K
     - Fe(C\ :sub:`5`\ H\ :sub:`5`\ )\ :sub:`2` (molecular)

Each material is tested in both **SCF** (self-consistent field) and **noSCF**
modes, for a total of 16 test cases.

Running the Benchmarks
----------------------

After building the project (see :doc:`installation`):

.. code-block:: bash

   python tests/integration/benchmark.py

The script will:

1. Run ``feff8l`` in each material's ``SCF/`` and ``noSCF/`` folders
2. Parse the resulting ``chi.dat`` files
3. Compare against the Fortran baseline in the ``baseline/`` folders
4. Compute R-factors for each comparison
5. Generate summary plots
6. Print a pass/fail summary table

Output Files
------------

The benchmark produces two summary plots in ``tests/integration/``:

- ``ALL_SCF_final.png`` --- Overlay of C++ vs Fortran :math:`\chi(k)` for all
  materials in SCF mode
- ``ALL_noSCF_final.png`` --- Same for noSCF mode

Each subplot shows the Fortran baseline and C++ result overlaid, with the
R-factor displayed in the title.

R-Factor Metric
---------------

The R-factor quantifies the agreement between the C++ result and the Fortran
baseline:

.. math::

   R = \frac{\sum_i [\chi_\mathrm{C++}(k_i) - \chi_\mathrm{Fortran}(k_i)]^2}
            {\sum_i [\chi_\mathrm{Fortran}(k_i)]^2}

A result is marked **PASS** if :math:`R < 0.01` (less than 1% relative
difference) and **CHECK** otherwise.

Test Directory Structure
------------------------

.. code-block:: text

   tests/integration/
   +-- benchmark.py          # Benchmark runner script
   +-- work/
       +-- Copper/
       |   +-- SCF/           # Contains feff.inp for SCF test
       |   +-- noSCF/         # Contains feff.inp for noSCF test
       |   +-- baseline/
       |       +-- withSCF/   # Fortran reference output (SCF)
       |       +-- noSCF/     # Fortran reference output (noSCF)
       +-- NiO/
       |   +-- ...
       +-- (6 more materials)

The ``SCF/`` and ``noSCF/`` folders contain only ``feff.inp`` files and are
populated with output when the benchmark runs. The ``baseline/`` folders
contain the reference Fortran output and should not be modified.

Prerequisites
-------------

The benchmark script requires Python 3 with:

.. code-block:: bash

   pip install numpy matplotlib
