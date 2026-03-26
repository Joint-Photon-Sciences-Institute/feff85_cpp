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

chi.dat --- EXAFS :math:`\chi(k)`
---------------------------------

The primary EXAFS output. Contains the sum of all significant scattering path
contributions. Comment lines begin with ``#``.

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Column
     - Units
     - Description
   * - ``k``
     - 1/Angstrom
     - Photoelectron wave vector
   * - ``chi``
     - (dimensionless)
     - :math:`\chi(k)` --- the EXAFS oscillation, computed as
       :math:`\mathrm{Im}\left[\sum_\Gamma \tilde{\chi}_\Gamma(k)\right]`
       where the sum runs over all included scattering paths
   * - ``mag``
     - (dimensionless)
     - :math:`|\chi(k)|` --- magnitude of the complex :math:`\chi(k)`
   * - ``phase``
     - radians
     - Phase of the complex :math:`\chi(k)`, with :math:`2\pi` jumps removed
       (unwrapped)

The ``chi`` column is the quantity most commonly used in EXAFS data analysis.
It is **not** :math:`k`-weighted; to plot :math:`k^n\chi(k)`, multiply the
``chi`` column by ``k^n``.

xmu.dat --- X-ray Absorption :math:`\mu(E)`
--------------------------------------------

The full X-ray absorption spectrum. Contains the total absorption coefficient
and its components, normalized to the cross-section at 50 eV above the edge.

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Column
     - Units
     - Description
   * - ``omega``
     - eV
     - Absolute photon energy
   * - ``e``
     - eV
     - Energy relative to the absorption edge (:math:`E - E_0`)
   * - ``k``
     - 1/Angstrom
     - Photoelectron wave vector :math:`k = \sqrt{2m(E-E_0)/\hbar^2}`
   * - ``mu``
     - (normalized)
     - :math:`\mu(E) = (\mu_0(E)(1 + \chi(E)) / \mu_0(E_0 + 50\,\mathrm{eV})`
       --- total absorption including fine structure, normalized to the
       atomic background cross-section evaluated at 50 eV above the edge
   * - ``mu0``
     - (normalized)
     - :math:`\mu_0(E) / \mu_0(E_0 + 50\,\mathrm{eV})` --- smooth atomic
       background absorption (no fine structure), same normalization
   * - ``chi``
     - (dimensionless)
     - :math:`\chi(E)` --- the fine-structure oscillation at this energy,
       equivalent to :math:`(\mu - \mu_0) / \mu_0`

The normalization reference :math:`\mu_0(E_0 + 50\,\mathrm{eV})` is the
atomic absorption cross-section interpolated at 50 eV above the edge energy.
This value is printed in the file header. This normalization ensures that
the edge step is approximately 1.0, making ``mu`` directly comparable across
different calculations.

feffNNNN.dat --- Individual Scattering Path Data
-------------------------------------------------

One file per significant scattering path (e.g., ``feff0001.dat``,
``feff0002.dat``, ...). These files contain the effective scattering amplitude
and phase for a single path, which can be used directly in the EXAFS equation
to compute that path's contribution to :math:`\chi(k)`.

**Header section:**

Each file begins with a header describing the path geometry:

- **nleg** --- number of legs in the scattering path
- **deg** --- path degeneracy :math:`N_\Gamma`
- **reff** --- effective half-path length :math:`R_\Gamma` (Angstroms)
- **rnrmav** --- average Norman radius (Bohr)
- **edge** --- absorption edge energy (eV)
- Atom positions along the path (Angstroms), with potential type and
  atomic number

**Data columns:**

.. list-table::
   :header-rows: 1
   :widths: 18 15 67

   * - Column
     - Units
     - Description
   * - ``k``
     - 1/Angstrom
     - Photoelectron wave vector
   * - ``real[2*phc]``
     - radians
     - :math:`2\delta_c + l\pi` --- twice the central-atom phase shift plus
       :math:`l\pi` (where :math:`l` is the angular momentum of the final
       state). This is the phase contribution from the absorbing atom.
   * - ``mag[feff]``
     - Angstroms
     - :math:`|f_\mathrm{eff}(k)|` --- magnitude of the effective
       curved-wave scattering amplitude. This is the amplitude that appears
       in the EXAFS equation **after** removing the :math:`1/kR^2` geometric
       factor, the mean free path, and the reduction factor.
   * - ``phase[feff]``
     - radians
     - :math:`\arg[f_\mathrm{eff}(k)]` --- phase of the effective scattering
       amplitude (with the central-atom phase already removed)
   * - ``red factor``
     - (dimensionless)
     - :math:`\exp(-2\,\mathrm{Im}[\delta_c])` --- reduction factor from
       inelastic losses at the central atom
   * - ``lambda``
     - Angstroms
     - :math:`\lambda(k) = 1/\mathrm{Im}[k]` --- photoelectron mean free path
   * - ``real[p]``
     - 1/Angstrom
     - :math:`\mathrm{Re}[p(k)]` --- real part of the complex momentum
       (photoelectron wavenumber including the real self-energy shift)

**Reconstructing** :math:`\chi(k)` **from feffNNNN.dat:**

The contribution of a single path to :math:`\chi(k)` can be computed from the
columns as:

.. math::

   \chi_\Gamma(k) = S_0^2 \cdot \frac{N_\Gamma \cdot |f_\mathrm{eff}(k)|}{k \cdot R_\Gamma^2}
   \cdot \sin\!\Big(2kR_\Gamma + \mathrm{real}[2\varphi_c] + \mathrm{phase}[f_\mathrm{eff}]\Big)
   \cdot e^{-2\sigma_\Gamma^2 k^2}
   \cdot e^{-2R_\Gamma/\lambda(k)}

where :math:`N_\Gamma` and :math:`R_\Gamma` come from the header, ``mag[feff]``
provides :math:`|f_\mathrm{eff}|`, ``phase[feff]`` provides
:math:`\arg[f_\mathrm{eff}]`, ``real[2*phc]`` provides
:math:`2\delta_c + l\pi`, ``lambda`` provides :math:`\lambda(k)`, and
:math:`\sigma^2` is the Debye-Waller factor (computed separately or supplied
via the ``DEBYE`` or ``SIG2`` cards).

These files are the standard input for EXAFS fitting programs such as
`Artemis/IFEFFIT <https://bruceravel.github.io/demeter/>`_ and
`Larch <https://xraypy.github.io/xraylarch/>`_, which refine structural
parameters (:math:`N`, :math:`R`, :math:`\sigma^2`, :math:`\Delta E_0`) by
fitting the theoretical :math:`\chi(k)` to experimental data.
