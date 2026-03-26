#pragma once
// Type definitions for ATOM module — replaces Fortran COMMON blocks.
// These structs hold the shared mutable state passed between ATOM subroutines.
//
// ATOM works with real arrays on a 251-point radial grid.
// FOVRG uses complex versions on the larger nrptx-point grid.

#include <feff/types.hpp>
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cstring>

namespace feff::atom {

/// Maximum number of orbitals
inline constexpr int max_orb = 30;

/// ATOM radial grid size (must be odd for Simpson integration)
inline constexpr int atom_grid = 251;

/// Maximum development coefficients at origin
inline constexpr int max_dev = 10;

// =========================================================================
// Orbital wave functions and development coefficients
// Replaces unnamed COMMON (cg,cp,bg,bp,fl,fix,ibgp) in ATOM
// and COMMON /dff/ in FOVRG
// =========================================================================

/// Real orbital arrays for ATOM (251-point grid)
struct OrbitalArraysReal {
    double cg[atom_grid][max_orb] = {};   // large radial component
    double cp[atom_grid][max_orb] = {};   // small radial component
    double bg[max_dev][max_orb] = {};     // large component dev. coefficients
    double bp[max_dev][max_orb] = {};     // small component dev. coefficients
    double fl[max_orb] = {};              // power of origin behavior
    double fix[max_orb] = {};             // quick fix factors
    int ibgp = 0;
};

/// Real orbital arrays for FOVRG (nrptx-point grid)
struct OrbitalArraysFovrg {
    double cg[nrptx][max_orb] = {};
    double cp[nrptx][max_orb] = {};
    double bg[max_dev][max_orb] = {};
    double bp[max_dev][max_orb] = {};
    double fl[max_orb] = {};
    double fix[max_orb] = {};
    int ibgp = 0;
};

// =========================================================================
// Orbital configuration — shared between ATOM and FOVRG
// Replaces COMMON /ratom1/
// =========================================================================
struct OrbitalConfig {
    double xnel[max_orb] = {};    // occupation numbers
    double en[max_orb] = {};      // orbital energies (Hartrees)
    double scc[max_orb] = {};     // convergence acceleration factors
    double scw[max_orb] = {};     // wave function convergence errors
    double sce[max_orb] = {};     // energy convergence errors
    int nq[max_orb] = {};         // principal quantum numbers
    int kap[max_orb] = {};        // kappa quantum numbers
    int nmax[max_orb] = {};       // last tabulation point for each orbital
};

// =========================================================================
// SCF iteration parameters
// Replaces COMMON /itescf/
// =========================================================================
struct ScfParams {
    double testy = 1.0e-5;       // wave function precision criterion
    double rap[2] = {};           // convergence ratios
    double teste = 5.0e-6;       // energy precision criterion
    int nz = 0;                   // nuclear charge Z
    int norb = 0;                 // total number of orbitals
    int norbsc = 0;               // number of core orbitals
};

// =========================================================================
// Dirac equation working arrays — real version for ATOM (251 pts)
// Replaces COMMON /comdir/
// =========================================================================
struct DiracWorkspaceReal {
    double cl = 0.0;             // speed of light (a.u.)
    double dz = 0.0;             // nuclear charge
    double dg[atom_grid] = {};   // large component workspace
    double ag[max_dev] = {};     // development coefficients for dg
    double dp[atom_grid] = {};   // small component workspace
    double ap[max_dev] = {};     // development coefficients for dp
    double eg[atom_grid] = {};   // exchange potential (large component)
    double ceg[max_dev] = {};    // dev. coefficients for eg
    double ep[atom_grid] = {};   // exchange potential (small component)
    double cep[max_dev] = {};    // dev. coefficients for ep
    double dv[atom_grid] = {};   // direct potential (Coulomb + nuclear)
    double av[max_dev] = {};     // dev. coefficients for dv
};

// =========================================================================
// Dirac equation working arrays — complex version for FOVRG (nrptx pts)
// Replaces COMMON /comdic/
// =========================================================================
struct DiracWorkspaceComplex {
    double cl = 0.0;             // speed of light (a.u.)
    double dz = 0.0;             // nuclear charge
    FeffComplex gg[nrptx] = {};  // large component (regular solution)
    FeffComplex ag[max_dev] = {};
    FeffComplex gp[nrptx] = {};  // small component (regular solution)
    FeffComplex ap[max_dev] = {};
    FeffComplex dv[nrptx] = {};  // direct potential
    FeffComplex av[max_dev] = {};
    FeffComplex eg[nrptx] = {};  // exchange potential (large)
    FeffComplex ceg[max_dev] = {};
    FeffComplex ep[nrptx] = {};  // exchange potential (small)
    FeffComplex cep[max_dev] = {};
};

// =========================================================================
// Lagrange parameters and shell data
// Replaces COMMON /scrhf1/
// =========================================================================
struct LagrangeParams {
    // Lagrange multipliers stored with triangular indexing:
    // eps[i + (j-1)*(j-2)/2] for orbital pair (i,j)
    double eps[435] = {};
    int nre[max_orb] = {};       // shell occupation flags
    int ipl = 0;                 // open-shell flag
};

// =========================================================================
// Nuclear potential and development
// Real version for ATOM (251 pts), replaces COMMON /snoyau/
// =========================================================================
struct NuclearPotentialReal {
    double dvn[atom_grid] = {};  // nuclear potential
    double anoy[max_dev] = {};   // development coefficients at origin
    int nuc = 0;                 // nuclear radius index (1=point charge)
};

// =========================================================================
// Nuclear potential — complex version for FOVRG (nrptx pts)
// Replaces COMMON /snoyac/
// =========================================================================
struct NuclearPotentialComplex {
    double dvn[nrptx] = {};      // nuclear potential (still real)
    double anoy[max_dev] = {};
    int nuc = 0;
};

// =========================================================================
// Radial mesh and test parameters
// Real version for ATOM (251 pts), replaces COMMON /tabtes/
// =========================================================================
struct MeshParamsReal {
    double hx = 0.0;            // exponential mesh step dx
    double dr[atom_grid] = {};   // radial mesh points r(i) = exp((i-1)*hx - 8.8)
    double test1 = 0.0;
    double test2 = 0.0;
    int ndor = 10;               // number of development coefficients
    int np = 0;                  // number of tabulation points
    int nes = 0;
    int method = 0;              // integration method (1 or 2)
    int idim = atom_grid;        // grid dimension
};

// =========================================================================
// Radial mesh — complex version for FOVRG (nrptx pts)
// Replaces COMMON /tabtec/
// =========================================================================
struct MeshParamsComplex {
    double hx = 0.0;
    double dr[nrptx] = {};
    double test1 = 0.0;
    double test2 = 0.0;
    int ndor = 10;
    int np = 0;
    int nes = 0;
    int method = 0;
    int idim = nrptx;
};

// =========================================================================
// Angular Coulomb coefficients — real for ATOM
// Replaces COMMON /mulabk/
// =========================================================================
struct AngularCoefficients {
    // afgk(i,j,k) for orbitals i,j and multipole k=0..3
    double afgk[max_orb][max_orb][4] = {};
};

// =========================================================================
// Angular exchange coefficients — for FOVRG
// Replaces COMMON /mulabc/
// afgkc(-ltot-1:ltot, 30, 0:3) → C++ offset: index + ltot + 1
// =========================================================================
struct AngularCoefficientsC {
    static constexpr int offset = ltot + 1;
    static constexpr int dim1 = 2 * ltot + 2;  // -ltot-1..ltot
    double afgkc[dim1][max_orb][4] = {};

    double& operator()(int i, int j, int k) {
        return afgkc[i + offset][j][k];
    }
    double operator()(int i, int j, int k) const {
        return afgkc[i + offset][j][k];
    }
};

// =========================================================================
// Breit term coefficients (ATOM only)
// Replaces COMMON /tabre/
// =========================================================================
struct BreitCoefficients {
    double cmag[3] = {};         // magnetic coefficients
    double cret[3] = {};         // retardation coefficients
};

// =========================================================================
// Dirac solver internal state (ATOM only)
// Replaces COMMON /subdir/
// =========================================================================
struct DiracSolverState {
    double ell = 0.0;
    double fk = 0.0;
    double ccl = 0.0;
    int imm = 0;
    int nd = 0;
    int node = 0;
    int mat = 0;
};

// =========================================================================
// Error state
// Replaces COMMON /messag/
// =========================================================================
struct ErrorState {
    char dlabpr[9] = {};         // label for error message
    int numerr = 0;              // error number
};

// =========================================================================
// Inelastic flag
// Replaces COMMON /inelma/
// =========================================================================
struct InelasticFlag {
    int nem = 0;                 // 0=full density, !=0=exchange integrals
};

// =========================================================================
// Print control
// Replaces COMMON /print/
// =========================================================================
struct PrintControl {
    int iprint = 0;
};

// =========================================================================
// Aggregate state for ATOM SCF calculation
// =========================================================================
struct AtomState {
    OrbitalArraysReal orb;       // cg, cp, bg, bp, fl, fix
    OrbitalConfig config;        // xnel, en, nq, kap, etc.
    ScfParams scf;               // iteration control
    DiracWorkspaceReal work;     // working arrays
    LagrangeParams lagrange;     // eps, nre, ipl
    NuclearPotentialReal nuclear; // dvn, anoy, nuc
    MeshParamsReal mesh;         // hx, dr, etc.
    AngularCoefficients angular; // afgk
    BreitCoefficients breit;     // cmag, cret
    DiracSolverState solver;     // ell, fk, etc.
    ErrorState error;            // dlabpr, numerr
    InelasticFlag inelastic;     // nem
    PrintControl print;          // iprint
};

// =========================================================================
// Aggregate state for FOVRG calculation
// =========================================================================
struct FovrgState {
    OrbitalArraysFovrg orb;      // cg, cp, bg, bp, fl, fix (nrptx grid)
    OrbitalConfig config;        // shared format with ATOM
    ScfParams scf;
    DiracWorkspaceComplex work;  // complex working arrays
    LagrangeParams lagrange;
    NuclearPotentialComplex nuclear;
    MeshParamsComplex mesh;
    AngularCoefficientsC angular; // afgkc with offset indexing
    ErrorState error;
};

} // namespace feff::atom
