#pragma once
// Complex radial integral functions for FOVRG module.
// Converted from: aprdec.f, aprdep.f, yzktec.f, yzkrdc.f, dsordc.f, diff.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;
using feff::atom::OrbitalArraysFovrg;
using feff::atom::MeshParamsComplex;
using feff::atom::DiracWorkspaceComplex;

/// Complex polynomial product: coefficient of power (l-1) for two polynomials.
/// @param a  Complex polynomial coefficients (size 10)
/// @param b  Real polynomial coefficients (size 10)
/// @param l  Power index (1-based term)
/// @return   Sum of a[m]*b[l-m] for m=1..l
FeffComplex aprdec(const FeffComplex a[10], const double b[10], int l);

/// Real polynomial product: coefficient of power (l-1) for two polynomials.
/// @param a  Real polynomial coefficients (size 10)
/// @param b  Real polynomial coefficients (size 10)
/// @param l  Power index (1-based term)
/// @return   Sum of a[m]*b[l-m] for m=1..l
double aprdep(const double a[10], const double b[10], int l);

/// Calculate yk(r) = zk(r) + r^(k+1) * integral(r..inf, f(u)*u^(-k-1)).
/// Integration from point to point by a 4-point method.
///
/// @param f     Input function / output yk (size nrptx)
/// @param af    Development coefficients of f / output yk coefficients (size 10)
/// @param g     Output zk (size nrptx)
/// @param ag    Output zk development coefficients (size 10)
/// @param dr    Radial mesh (size nrptx)
/// @param ap    Power of first term / output constant for yk [in/out]
/// @param h     Exponential step
/// @param k     Multipole order
/// @param nd    Number of development terms
/// @param np    Number of tabulation points [modified: min(np, idim-1)]
/// @param idim  Array dimension
/// @param dyzk  Boundary term for yk at np+1
void yzktec(FeffComplex f[], FeffComplex af[10], FeffComplex g[],
            FeffComplex ag[10], const double dr[], FeffComplex& ap,
            double h, int k, int nd, int& np, int idim, FeffComplex& dyzk);

/// Calculate function yk using orbital products.
/// yk = r * integral of f(s)*uk(r,s) where uk is the Green function kernel.
///
/// @param i     Orbital index (0-based)
/// @param k     Multipole order
/// @param flps  Power of photoelectron wave function
/// @param ps    Photoelectron large component (size nrptx)
/// @param qs    Photoelectron small component (size nrptx)
/// @param aps   Development coefficients for ps (size 10)
/// @param aqs   Development coefficients for qs (size 10)
/// @param state FOVRG state containing orbitals, mesh, workspace
void yzkrdc(int i, int k, double flps, const FeffComplex ps[],
            const FeffComplex qs[], const FeffComplex aps[10],
            const FeffComplex aqs[10], FovrgState& state);

/// Calculate overlap integral by Simpson method.
/// Integrates hg(r)*r where hg(l) = dg(l)*cg(l,j) + dp(l)*cp(l,j).
///
/// @param j     Orbital index (0-based)
/// @param a     Power of first term of dg, dp at origin
/// @param dg    Large component array (size nrptx)
/// @param dp    Small component array (size nrptx)
/// @param ag    Development coefficients of dg (size 10)
/// @param ap    Development coefficients of dp (size 10)
/// @param orb   Orbital arrays
/// @param mesh  Mesh parameters
/// @return      Overlap integral value
FeffComplex dsordc(int j, FeffComplex a, const FeffComplex dg[],
                   const FeffComplex dp[], const FeffComplex ag[10],
                   const FeffComplex ap[10], OrbitalArraysFovrg& orb,
                   MeshParamsComplex& mesh);

/// Calculate vm(i) = (dV/dx)*r(i)*(kap+1)/cl for the c3 correction term.
/// Reference: Koelling, Harmon, J.Phys.C, 3107 (1977), eq. 14.
///
/// @param v    Input potential (complex, size n)
/// @param dr   Radial mesh (size n)
/// @param kap  Kappa quantum number
/// @param cl   Speed of light
/// @param dx   Exponential step
/// @param n    Number of points
/// @param vm   Output array (size n)
void diff(const FeffComplex v[], const double dr[], int kap,
          double cl, double dx, int n, FeffComplex vm[]);

} // namespace feff::fovrg
