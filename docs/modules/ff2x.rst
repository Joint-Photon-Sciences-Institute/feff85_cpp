FF2X Module
===========

**Namespace:** ``feff::ff2x``

**Header:** ``src/ff2x/ff2x.hpp``

The FF2X module generates the final output spectra by summing contributions
from all scattering paths. This is Stage 6 of the pipeline.

Overview
--------

FF2X reads the effective scattering amplitudes from ``feff.pad`` and the
cross-section data from ``xsect.json``, then combines them with Debye-Waller
factors and the mean free path to produce the output spectra.

Three output modes are supported:

- **EXAFS** --- :math:`\chi(k)` via ``ff2chi``
- **XANES** --- :math:`\mu(E)` via ``ff2xmu``
- **Anomalous** --- anomalous scattering / DANES via ``ff2afs``

Key Functions
-------------

``void ff2chi(const FF2xParams& p, int iabs)``
   Calculate EXAFS :math:`\chi(k)` by summing path contributions. Applies
   Debye-Waller factors, mean free path, and :math:`S_0^2`.

``void ff2xmu(const FF2xParams& p, int iabs)``
   Calculate XANES :math:`\mu(E)` including the atomic background.

``void ff2afs(const FF2xParams& p, int iabs)``
   Calculate anomalous scattering factors.

``void xscorr(...)``
   Apply self-energy corrections to :math:`\chi(k)`.

``void dwadd(...)``
   Add Debye-Waller factors to path contributions using the correlated
   Debye model.

Data Structures
---------------

``FF2xParams``
   Input parameters including :math:`S_0^2`, Debye temperature, temperature,
   output type (EXAFS/XANES), and spectral broadening.

``FeffPadData``
   Data read from ``feff.pad``: path amplitudes, phases, and geometry.

``XsectData``
   Cross-section and energy mesh from ``xsect.json``.

``PathListEntry``
   Individual path from ``list.dat`` with index, degeneracy, and path length.

Output
------

- ``chi.dat`` --- EXAFS :math:`\chi(k)` (k, chi, magnitude, phase)
- ``xmu.dat`` --- X-ray absorption :math:`\mu(E)` vs energy
