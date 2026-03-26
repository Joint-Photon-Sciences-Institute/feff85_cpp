Math Module
===========

**Namespace:** ``feff::math``

**Headers:** ``src/math/``

The Math module provides mathematical functions used throughout the FEFF
calculation. All functions operate on double-precision values.

Bessel Functions
----------------

**Header:** ``src/math/bessel.hpp``

Spherical Bessel functions :math:`j_l(x)`, :math:`n_l(x)`, and spherical
Hankel functions :math:`h_l^{(1)}(x)` for complex arguments.

Wigner Symbols
--------------

**Header:** ``src/math/wigner.hpp``

Wigner 3j and 6j symbols, Clebsch-Gordan coefficients, and rotation matrix
elements (Wigner d-functions) used in angular momentum coupling.

Legendre Polynomials
--------------------

**Header:** ``src/math/legendre.hpp``

Legendre polynomials :math:`P_l(x)` and associated Legendre functions
:math:`P_l^m(x)` used in multipole expansions.

Interpolation
-------------

**Header:** ``src/math/interpolation.hpp``

Cubic spline interpolation and linear interpolation routines for
transferring data between grids.

Integration
-----------

**Header:** ``src/math/integration.hpp``

Numerical integration on the logarithmic radial grid using Simpson's rule
and related quadrature methods.

Other Functions
---------------

**Header:** ``src/math/determinant.hpp``
   Matrix determinants for complex matrices.

**Header:** ``src/math/polynomial_roots.hpp``
   Root-finding for polynomials.

**Header:** ``src/math/convolution.hpp``
   Convolution routines for spectral broadening.

**Header:** ``src/math/distance.hpp``
   Interatomic distance calculations.

**Header:** ``src/math/phase_amplitude.hpp``
   Phase and amplitude extraction from complex scattering amplitudes.

**Header:** ``src/math/sommerfeld.hpp``
   Sommerfeld expansion coefficients for thermal corrections.

**Header:** ``src/math/bcoef.hpp``
   Binomial and related combinatorial coefficients.
