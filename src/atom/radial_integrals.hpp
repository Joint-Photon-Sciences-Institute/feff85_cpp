#pragma once
// Radial integral routines for the ATOM module.
// Converted from: yzkteg.f, yzkrdf.f, dsordf.f, fdrirk.f, fdmocc.f, akeato.f

#include "atom_types.hpp"

namespace feff::atom {

/// Core yk/zk integration kernel.
/// Computes yk(r) = zk(r) + r^(k+1) * integral_{r}^{inf} f(u)*u^(-k-1) du
/// where  zk(r) = r^(-k) * integral_{0}^{r} f(u)*u^k du.
///
/// On entry: f[] contains the input function, af[] its dev. coefficients,
///           ap is the power of origin behavior.
/// On exit:  f[] contains yk, g[] contains zk, af[]/ag[] hold dev. coefficients,
///           ap holds the constant for yk at origin.
///
/// All arrays are 0-based. h = exponential mesh step, k = multipole order,
/// nd = number of dev. coefficients, np = number of tabulation points,
/// idim = array dimension.
void yzkteg(double f[], double af[], double g[], double ag[], const double dr[],
            double& ap, double h, int k, int nd, int np, int idim);

/// Calculate yk function for orbital products.
/// i, j = orbital indices (1-based Fortran convention, or i<=0 for pre-constructed f).
/// When i<=0, j is a grid point count (not an orbital index).
/// k = multipole order.
/// Results stored in work.dp (yk) and work.dg (zk).
void yzkrdf(int i, int j, int k,
            OrbitalArraysReal& orb, DiracWorkspaceReal& work,
            OrbitalConfig& config, InelasticFlag& inel, MeshParamsReal& mesh);

/// Simpson overlap integral of hg(r)*r^n.
/// i, j = orbital indices (0-based), n = power of r, jnd = construction mode.
/// a = power of origin behavior (for jnd >= 4 or negative jnd).
/// Returns the integral value.
double dsordf(int i, int j, int n, int jnd, double a,
              OrbitalArraysReal& orb, DiracWorkspaceReal& work,
              OrbitalConfig& config, MeshParamsReal& mesh);

/// Radial integrals Rk = integral f(r)*uk(r,s)*g(s).
/// i,j define f; l,m define g; k = multipole order.
/// All orbital indices are 0-based.
double fdrirk(int i, int j, int l, int m, int k,
              OrbitalArraysReal& orb, DiracWorkspaceReal& work,
              OrbitalConfig& config, InelasticFlag& inel, MeshParamsReal& mesh);

/// Product of occupation numbers for orbitals i and j (0-based).
/// Accounts for self-interaction correction when i == j.
double fdmocc(int i, int j, OrbitalConfig& config);

/// Direct Coulomb angular coefficient a_k(i,j).
/// i, j = orbital indices (0-based), k = multipole order (even).
/// Reads from ang.afgk[min(i,j)][max(i,j)][k/2].
double akeato(int i, int j, int k, AngularCoefficients& ang);

/// Exchange Coulomb angular coefficient b_k(i,j).
/// i, j = orbital indices (0-based), k = multipole order (even).
/// Reads from ang.afgk[max(i,j)][min(i,j)][k/2].
/// Returns 0 when i == j.
double bkeato(int i, int j, int k, AngularCoefficients& ang);

} // namespace feff::atom
