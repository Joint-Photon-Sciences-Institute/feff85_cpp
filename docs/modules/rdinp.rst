RDINP Module
============

**Namespace:** ``feff::rdinp``

**Header:** ``src/rdinp/rdinp.hpp``

The RDINP module parses the ``feff.inp`` input file and populates the
``FeffInput`` structure. This is the first stage of the pipeline.

Functions
---------

``int rdinp(FeffInput& inp, const std::string& file)``
   Parse the specified input file and populate ``inp`` with all calculation
   parameters. Returns 0 on success.

FeffInput Structure
-------------------

Defined in ``include/feff/feff_input.hpp``, the ``FeffInput`` struct is the
central data container holding all parameters from ``feff.inp``. It is
organized into sections corresponding to each pipeline module:

**Version Information:**
   ``vfeff``, ``vf85e`` --- version strings

**Atomic Structure (ATOMS card):**
   - ``ratx[natx][3]`` --- atomic coordinates in Angstroms
   - ``iphatx[natx]`` --- potential type index for each atom
   - ``nat`` --- number of atoms

**Global Parameters:**
   - ``polarization`` --- polarization vector (if enabled)
   - ``ellipticity`` --- ellipticity parameters
   - ``ispin``, ``nspx`` --- spin configuration

**MOD1 (POT) Parameters:**
   - ``nph`` --- number of unique potentials
   - ``iz[nphx+1]`` --- atomic numbers
   - ``nohole`` --- core-hole treatment flag
   - ``rfms1`` --- FMS radius for SCF
   - ``lfms1`` --- max angular momentum for FMS in SCF

**MOD2 (XSPH) Parameters:**
   - ``ixc``, ``vr0``, ``vi0`` --- exchange model and corrections
   - ``rgrd`` --- radial grid ratio

**MOD3 (FMS) Parameters:**
   - ``rfms2`` --- FMS radius
   - ``lfms2`` --- max angular momentum for FMS

**MOD4 (PATH) Parameters:**
   - ``rmax`` --- max half-path length
   - ``nleg`` --- max legs per path
   - ``pcritk``, ``pcrith`` --- path filtering criteria

**MOD5 (GENFMT) Parameters:**
   - ``iorder`` --- approximation order (2 = standard)
   - ``critcw`` --- curved-wave filter percentage

**MOD6 (FF2X) Parameters:**
   - ``s02`` --- amplitude reduction factor :math:`S_0^2`
   - ``tk``, ``thetad`` --- temperature and Debye temperature
   - ``ispec`` --- output type (EXAFS, XANES, etc.)
