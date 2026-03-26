Constants & Dimensions
======================

Physical Constants
------------------

Defined in ``include/feff/constants.hpp``:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Name
     - Value
     - Description
   * - ``pi``
     - 3.14159265...
     - :math:`\pi`
   * - ``coni``
     - (0, 1)
     - Complex unit :math:`i`
   * - ``bohr``
     - 0.52917721067
     - Bohr radius in Angstroms
   * - ``ryd``
     - 13.60569301
     - Rydberg energy in eV
   * - ``hart``
     - 27.21138602
     - Hartree energy in eV (= 2 Ryd)
   * - ``alpinv``
     - 137.035999139
     - Inverse fine-structure constant :math:`1/\alpha`
   * - ``alphfs``
     - 0.00729735...
     - Fine-structure constant :math:`\alpha`
   * - ``fa``
     - 1.91916...
     - Fermi wavevector prefactor :math:`k_F = f_a / r_s`
   * - ``raddeg``
     - 57.29578...
     - Radians to degrees conversion :math:`180/\pi`

Type Aliases
------------

Defined in ``include/feff/types.hpp``:

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Alias
     - Type
     - Description
   * - ``FeffComplex``
     - ``std::complex<double>``
     - Complex number type used throughout
   * - ``Vec3``
     - ``std::array<double, 3>``
     - 3D vector

Dimension Parameters
--------------------

Defined in ``include/feff/dimensions.hpp``. These control the maximum sizes of
internal arrays and must be set at compile time.

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   * - Name
     - Value
     - Description
   * - ``nclusx``
     - 100
     - Max atoms for FMS cluster
   * - ``natx``
     - 1000
     - Max atoms for path finder
   * - ``nphx``
     - 11
     - Max unique potential types (must be odd)
   * - ``ltot``
     - 24
     - Max angular momentum quantum number
   * - ``lx``
     - 4
     - Max orbital momentum for FMS
   * - ``nrptx``
     - 1251
     - Loucks radial grid size
   * - ``nex``
     - 150
     - Max energy points for EXAFS
   * - ``lamtot``
     - 15
     - Max distinct :math:`\lambda` values for GENFMT
   * - ``mtot``
     - 4
     - Max :math:`m` quantum number
   * - ``ntot``
     - 2
     - Max :math:`n` quantum number
   * - ``npatx``
     - 8
     - Max atoms in a scattering path
   * - ``legtot``
     - 9
     - Max legs in a scattering path (npatx + 1)
   * - ``novrx``
     - 8
     - Max overlap shells
   * - ``nheadx``
     - 30
     - Max header lines in input
   * - ``nspx``
     - 1
     - Max spins (1 = spin-averaged)
   * - ``MxPole``
     - 1000
     - Max poles for HL multipole self-energy
