JSON I/O Module
===============

**Namespace:** ``feff::json_io``

**Header:** ``include/feff/json_io.hpp``

The JSON I/O module handles serialization and deserialization of calculation
parameters between pipeline stages. It uses the
`nlohmann/json <https://github.com/nlohmann/json>`_ library.

Overview
--------

Each pipeline stage writes its configuration to JSON files after parsing, and
subsequent stages read the JSON files they need. This provides a clean
interface between modules and enables restarting calculations from intermediate
stages.

Writer Functions
----------------

Called by RDINP after parsing ``feff.inp``:

``void write_atoms_json(const FeffInput& inp)``
   Write atomic structure data (coordinates, potential types).

``void write_global_json(const FeffInput& inp, int nabs)``
   Write global parameters (polarization, spin, title lines).

``void write_geom_json(const FeffInput& inp, int nat, const double rat[][3], const int iphat[], const int iatph[])``
   Write geometry data (transformed coordinates in Bohr).

``void write_pot_json(const FeffInput& inp)``
   Write POT module parameters (SCF, overlap, exchange).

``void write_xsph_json(const FeffInput& inp)``
   Write XSPH module parameters (phase shift options).

``void write_path_json(const FeffInput& inp)``
   Write PATH module parameters (rmax, nleg, criteria).

``void write_genfmt_json(const FeffInput& inp)``
   Write GENFMT module parameters (iorder, critcw).

``void write_ff2x_json(const FeffInput& inp)``
   Write FF2X module parameters (S02, Debye, output type).

``void write_libpotph_json(const FeffInput& inp)``
   Write combined potentials and phase shift input for the libpotph library.

Called by computational stages:

``void write_xsect_json(...)``
   Write cross-section and energy mesh data from XSPH.

``void write_feff_json(...)``
   Write scattering path data from GENFMT (one per path).

Reader Functions
----------------

``void read_atoms_json(int& nat, double rat[][3], int iphat[])``
   Read atomic structure from ``atoms.json``.

``void read_global_json(...)``
   Read global parameters from ``global.json``.

``void read_geom_json(int& nat, int& nph, int iatph[], double rat[][3], int iphat[], int ibounc[])``
   Read geometry from ``geom.json``.

``void read_pot_json(...)``
   Read all POT module parameters from ``pot.json``.

``void read_xsect_json(int& ntit, std::string titles[], double& s02, ...)``
   Read cross-section data from ``xsect.json``.

``void read_titles_json(int& ntit, std::string titles[])``
   Read only the title lines from ``xsect.json``.

``void read_libpotph_json(...)``
   Read combined potentials and phases input from ``libpotph.json``.
