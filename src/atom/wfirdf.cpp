// Initial Thomas-Fermi orbital calculation via Dirac equation.
// Converted from: src/ATOM/wfirdf.f
//
// INDEX CONVENTION: Loop indices remain 1-based as in Fortran.
// Array accesses use the F() macro for 0-based C++ arrays.

#include "wfirdf.hpp"
#include "soldir.hpp"
#include "nucdev.hpp"
#include "utility.hpp"
#include "../common/logging.hpp"
#include <cmath>
#include <sstream>

namespace feff::atom {

// 1-based to 0-based index helper
#define F(arr, i) (arr)[(i) - 1]

void wfirdf(double en[], double ch, int nq[], int kap[], int nmax[],
            int ido, AtomState& state)
{
    // Aliases into AtomState
    auto& orb    = state.orb;
    auto& work   = state.work;
    auto& scf    = state.scf;
    auto& inelas = state.inelastic;
    auto& error  = state.error;
    auto& nuclear = state.nuclear;
    auto& mesh   = state.mesh;

    int norb = scf.norb;
    int nz   = scf.nz;

    // Speed of light in atomic units
    work.cl = 1.370373e+02;
    work.dz = nz;

    // Make r-mesh and calculate nuclear potential
    double hx = 5.0e-02;
    // Match Fortran: dr1 = nz * exp(-8.8), where -8.8 is single-precision
    // in Fortran. exp(single) returns single, then integer*single = single,
    // which is promoted to double on assignment. We replicate this exactly.
    double dr1 = static_cast<double>(nz * std::exp(-8.8f));
    mesh.hx = hx;

    nucdev(nuclear.anoy, mesh.dr, nuclear.dvn, work.dz, mesh.hx,
           nuclear.nuc, mesh.idim, mesh.ndor, dr1);

    // Compute fl and fix for each orbital
    double a = (work.dz / work.cl) * (work.dz / work.cl);
    if (nuclear.nuc > 1) a = 0.0;
    for (int j = 1; j <= norb; ++j) {
        double b = F(kap, j) * F(kap, j) - a;
        F(orb.fl, j) = std::sqrt(b);
        // Quick fix of development coefficients
        F(orb.fix, j) = std::pow(F(mesh.dr, 1), F(orb.fl, j) - std::abs(F(kap, j)));
    }

    // Calculate potential from Thomas-Fermi model
    for (int i = 1; i <= mesh.idim; ++i) {
        F(work.dv, i) = (dentfa(F(mesh.dr, i), work.dz, ch) + F(nuclear.dvn, i)) / work.cl;
    }
    if (error.numerr != 0) return;

    // Zero out exchange potentials
    for (int i = 1; i <= mesh.idim; ++i) {
        F(work.eg, i) = 0.0;
        F(work.ep, i) = 0.0;
    }
    for (int i = 1; i <= orb.ibgp; ++i) {
        F(work.ceg, i) = 0.0;
        F(work.cep, i) = 0.0;
        F(work.av, i) = F(nuclear.anoy, i) / work.cl;
    }
    // Add Thomas-Fermi contribution at nuclear radius
    F(work.av, 2) = F(work.av, 2) + dentfa(F(mesh.dr, nuclear.nuc), work.dz, ch) / work.cl;

    mesh.test1 = scf.testy / scf.rap[0];  // rap(1) in Fortran = rap[0] in C++
    double b_val = mesh.test1;

    // Only ido=1 is supported
    if (ido != 1) {
        feff::common::logger().wlog("only option ido=1 left");
        ido = 1;
    }

    // Solve Dirac equation for each orbital
    for (int j = 1; j <= norb; ++j) {
        // Set initial development coefficient for large component
        // bg(1,j) in Fortran -> orb.bg[j-1][0] in C++ (bg[max_dev][max_orb], stored as bg[dev][orb])
        // Wait -- in atom_types.hpp: bg[max_dev][max_orb], so bg[0][j-1] = bg(1,j)
        orb.bg[0][j - 1] = 1.0;
        int i_val = F(nq, j) - std::abs(F(kap, j));
        if (F(kap, j) < 0) i_val = i_val - 1;
        if (i_val % 2 == 0) orb.bg[0][j - 1] = -orb.bg[0][j - 1];

        if (F(kap, j) >= 0) {
            orb.bp[0][j - 1] = orb.bg[0][j - 1] * work.cl * (F(kap, j) + F(orb.fl, j)) / work.dz;
            if (nuclear.nuc > 1) orb.bg[0][j - 1] = 0.0;
        } else {
            orb.bp[0][j - 1] = orb.bg[0][j - 1] * work.dz / (work.cl * (F(kap, j) - F(orb.fl, j)));
            if (nuclear.nuc > 1) orb.bp[0][j - 1] = 0.0;
        }

        mesh.np = mesh.idim;
        // Fortran: en(j) = -dz*dz/nq(j)*nq(j)
        // Left-to-right evaluation: ((-dz*dz)/nq(j))*nq(j) = -dz^2
        // This is a known Fortran precedence quirk -- the intended hydrogen-like
        // energy -Z^2/n^2 is not what gets computed. We preserve the exact
        // Fortran behavior as a crude initial guess for soldir.
        F(en, j) = -work.dz * work.dz / F(nq, j) * F(nq, j);

        mesh.method = 0;
        int ifail = 0;
        double bg1j = orb.bg[0][j - 1];
        double bp1j = orb.bp[0][j - 1];

        // b_val (Fortran 'b') is passed as ainf to soldir. It persists across
        // orbital iterations — Fortran sets b=test1 ONCE before the loop and
        // lets soldir modify it through the ainf reference for each orbital.
        soldir(F(en, j), F(orb.fl, j), bg1j, bp1j, b_val,
               F(nq, j), F(kap, j), F(nmax, j), ifail,
               work, mesh, state.solver, error);

        if (error.numerr != 0) {
            messer(error);
            std::ostringstream slog;
            slog << "soldir failed in wfirdf for orbital nq,kappa "
                 << F(nq, j) << " " << F(kap, j);
            feff::common::logger().wlog(slog.str());
        } else {
            // Copy results back to orbital arrays
            // Fortran: bg(i,j) = ag(i), bp(i,j) = ap(i)
            for (int i = 1; i <= orb.ibgp; ++i) {
                orb.bg[i - 1][j - 1] = F(work.ag, i);
                orb.bp[i - 1][j - 1] = F(work.ap, i);
            }
            // Fortran: cg(i,j) = dg(i), cp(i,j) = dp(i)
            for (int i = 1; i <= mesh.np; ++i) {
                orb.cg[i - 1][j - 1] = F(work.dg, i);
                orb.cp[i - 1][j - 1] = F(work.dp, i);
            }
        }
    }

    inelas.nem = 0;
}

#undef F

} // namespace feff::atom
