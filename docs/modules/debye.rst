Debye Module
============

**Namespace:** ``feff::debye``

**Header:** ``src/debye/sigms.hpp``

The Debye module calculates Debye-Waller factors using the correlated Debye
model for thermal disorder in EXAFS.

Overview
--------

The Debye-Waller factor :math:`e^{-2\sigma^2 k^2}` accounts for thermal and
static disorder in the interatomic distances. FEFF uses the correlated Debye
model, which properly accounts for correlations in atomic motion along the
scattering path.

For a path :math:`\Gamma` with legs connecting atoms :math:`i_0, i_1, \ldots, i_N`,
the mean-square relative displacement is:

.. math::

   \sigma_\Gamma^2 = \sum_{j} \hat{R}_j \cdot \langle (u_{i_j} - u_{i_{j+1}})
   (u_{i_j} - u_{i_{j+1}}) \rangle \cdot \hat{R}_j

where :math:`u_i` are atomic displacement vectors.

Key Functions
-------------

``void sigms(...)``
   Calculate :math:`\sigma^2` for a single-scattering path using the
   correlated Debye model.

``void sigcl(...)``
   Calculate :math:`\sigma^2` for a general multiple-scattering path using
   the classical (high-temperature) limit.

``void sigem(...)``
   Einstein model for :math:`\sigma^2` (single-frequency approximation).

``void sigrm(...)``
   Raman model contribution to :math:`\sigma^2`.

Key Parameters
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Parameter
     - Description
   * - ``tk``
     - Temperature in Kelvin
   * - ``thetad``
     - Debye temperature in Kelvin
   * - ``sig2``
     - Calculated mean-square relative displacement (Angstroms\ :sup:`2`)
