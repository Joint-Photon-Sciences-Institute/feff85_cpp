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

``void write_atoms_json(...)``
   Write atomic structure data (coordinates, potential types).

``void write_global_json(...)``
   Write global parameters (polarization, spin, title lines).

``void write_geom_json(...)``
   Write geometry data (transformed coordinates in Bohr).

``void write_pot_json(...)``
   Write POT module parameters (SCF, overlap, exchange).

``void write_xsph_json(...)``
   Write XSPH module parameters (phase shift options).

``void write_path_json(...)``
   Write PATH module parameters (rmax, nleg, criteria).

``void write_genfmt_json(...)``
   Write GENFMT module parameters (iorder, critcw).

``void write_ff2x_json(...)``
   Write FF2X module parameters (S02, Debye, output type).

Called by computational stages:

``void write_xsect_json(...)``
   Write cross-section and energy mesh data from XSPH.

``void write_feff_json(...)``
   Write scattering path data from GENFMT.

Reader Functions
----------------

``FeffInput read_atoms_json(const std::string& path)``
   Read atomic structure from ``atoms.json``.

``void read_global_json(...)``
   Read global parameters from ``global.json``.

``void read_geom_json(...)``
   Read geometry from ``geom.json``.

``void read_pot_json(...)``
   Read all POT module parameters from ``pot.json``.

``XsectData read_xsect_json(...)``
   Read cross-section data from ``xsect.json``.
