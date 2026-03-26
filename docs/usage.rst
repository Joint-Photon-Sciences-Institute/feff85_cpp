Usage
=====

Running FEFF
------------

FEFF reads its input from a file named ``feff.inp`` in the current working
directory:

.. code-block:: bash

   cd /path/to/your/calculation
   feff8l

The program runs the complete 6-stage calculation pipeline and produces output
files in the same directory.

Command-Line Options
--------------------

.. code-block:: text

   feff8l [input_file] [options]

**Positional arguments:**

``input_file``
   Path to the input file. Defaults to ``feff.inp`` in the current directory.

**Options:**

``--skip-pot``
   Skip stages 1-2 (input parsing and potential calculation). Useful when
   re-running with existing ``pot.pad`` and ``phase.pad`` files.

``--skip-to-path``
   Skip stages 1-3 (input, potentials, and phase shifts). Start from path
   enumeration.

``--skip-to-genfmt``
   Skip stages 1-4. Start from F-matrix calculation using existing
   ``paths.dat``.

``--skip-to-ff2x``
   Skip stages 1-5. Only run the final spectra generation using existing
   ``feff.pad``.

Input File Format
-----------------

The ``feff.inp`` file uses a card-based format inherited from the original
Fortran FEFF code. Each card begins with a keyword and is followed by
parameters. Lines beginning with ``*`` are comments.

Key cards include:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Card
     - Description
   * - ``TITLE``
     - Title lines for the calculation
   * - ``EDGE``
     - Absorption edge (e.g., ``K``, ``L1``, ``L2``, ``L3``)
   * - ``POTENTIALS``
     - Defines unique potential types with atomic numbers
   * - ``ATOMS``
     - Atomic coordinates (x, y, z in Angstroms), potential type, and tag
   * - ``CONTROL``
     - Flags to enable/disable each pipeline stage
   * - ``SCF``
     - Self-consistent field iteration parameters
   * - ``RMAX``
     - Maximum half-path length for path enumeration (Angstroms)
   * - ``NLEG``
     - Maximum number of legs in a scattering path
   * - ``POLARIZATION``
     - Polarization vector for polarized EXAFS
   * - ``EXCHANGE``
     - Exchange-correlation potential model (0=HL, 1=DH, 2=GW, 5=DH+HL)
   * - ``S02``
     - Amplitude reduction factor :math:`S_0^2`
   * - ``DEBYE``
     - Temperature and Debye temperature for Debye-Waller factors
   * - ``ELLIPTICITY``
     - Ellipticity vector for elliptical polarization
   * - ``PRINT``
     - Verbosity control for each module

Output Files
------------

After a successful run, the following files are produced:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - Description
   * - ``chi.dat``
     - EXAFS :math:`\chi(k)` data (k, chi, magnitude, phase)
   * - ``xmu.dat``
     - X-ray absorption coefficient :math:`\mu(E)` vs energy
   * - ``feffNNNN.dat``
     - Individual scattering path contributions (one per significant path)
   * - ``paths.dat``
     - Summary of all scattering paths with degeneracies and distances
   * - ``files.dat``
     - List of output files and their descriptions
   * - ``feff.pad``
     - Binary pad file with combined path amplitudes and phases
   * - ``phase.pad``
     - Binary pad file with phase shift data
   * - ``pot.pad``
     - Binary pad file with potential data
   * - ``feff.log``
     - Execution log with timing information

Intermediate JSON files (``atoms.json``, ``global.json``, ``pot.json``, etc.)
are also written and used for inter-module communication.
