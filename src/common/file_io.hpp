#pragma once
// File I/O for PAD-format data files (pot.pad, phase.pad) and headers.
// Converted from: src/COMMON/rdpot.f, src/COMMON/rdxsph.f,
//                 src/COMMON/head.f, src/COMMON/rdhead.f, src/COMMON/rdcmt.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <array>

namespace feff::common {

// ---------------------------------------------------------------------------
// PotData -- everything read from pot.pad (replaces rdpot outputs)
// ---------------------------------------------------------------------------
struct PotData {
    // Title lines
    int ntitle = 0;
    std::array<std::string, nheadx> title{};

    // Scalar parameters (from the 13-element dum block)
    double rnrmav = 0.0;    // average Norman radius
    double xmu    = 0.0;    // Fermi level
    double vint   = 0.0;    // muffin-tin zero (interstitial potential)
    double rhoint = 0.0;    // interstitial density
    double emu    = 0.0;    // edge position
    double s02    = 0.0;    // many-body reduction factor
    double erelax = 0.0;    // relaxation energy
    double wp     = 0.0;    // plasmon frequency estimate
    double ecv    = 0.0;    // core-valence separation energy
    double rs     = 0.0;    // r_s from rhoint
    double xf     = 0.0;    // kf from rhoint
    double qtotel = 0.0;    // total cluster charge
    double totvol = 0.0;    // total volume

    // Integer header fields
    int nph    = 0;
    int npadx  = 0;
    int nohole = 0;
    int ihole  = 0;
    int inters = 0;
    int iafolp = 0;
    int jumprm = 0;
    int iunf   = 0;

    // Per-potential-type integer arrays  (0:nphx)
    std::array<int, nphx + 1>  imt{};
    std::array<int, nphx + 1>  inrm{};
    std::array<int, nphx + 1>  iz{};

    // Per-potential-type double arrays (0:nphx)
    std::array<double, nphx + 1> rmt{};
    std::array<double, nphx + 1> rnrm{};
    std::array<double, nphx + 1> folp{};
    std::array<double, nphx + 1> folpx{};
    std::array<double, nphx + 1> xnatph{};
    std::array<double, nphx + 1> xion{};
    std::array<double, nphx + 1> qnrm{};

    // Kappa array (1:30)
    std::array<int, 30> kappa{};

    // Dirac spinor components  (flat storage, Fortran column-major order)
    std::array<double, 251>                     dgc0{};   // (251)
    std::array<double, 251>                     dpc0{};   // (251)
    std::vector<double> dgc;    // 251*30*(nph+1)
    std::vector<double> dpc;    // 251*30*(nph+1)
    std::vector<double> adgc;   // 10*30*(nph+1)
    std::vector<double> adpc;   // 10*30*(nph+1)

    // Density / potential arrays   251*(nph+1)
    std::vector<double> edens;
    std::vector<double> vclap;
    std::vector<double> vtot;
    std::vector<double> edenvl;
    std::vector<double> vvalgs;
    std::vector<double> dmag;

    // Valence info  30*(nph+1)
    std::vector<double> xnval;

    // Orbital energies (1:30)
    std::array<double, 30> eorb{};

    // Orbital indices  iorb(-4:3, 0:nphx) stored as flat [8][nphx+1]
    // Index mapping: iorb[i+4][iph]
    std::vector<int> iorb;   // 8*(nph+1)

    // SCF occupation numbers  xnmues(0:lx, 0:nphx) = (lx+1)*(nph+1)
    std::vector<double> xnmues;
};

// ---------------------------------------------------------------------------
// PhaseData -- everything read from phase.pad (replaces rdxsph outputs)
// ---------------------------------------------------------------------------
struct PhaseData {
    // Energy grid
    int ne  = 0;   // total energy points
    int ne1 = 0;   // points on main horizontal axis
    int ne3 = 0;   // points on auxiliary horizontal axis
    int nph = 0;   // number of potential types
    int ihole = 0;
    int ik0   = 0; // grid index at Fermi level
    int npadx = 0;
    int ixc   = 0; // potential model
    int nsp   = 0; // number of spins

    double rnrmav = 0.0;
    double xmu    = 0.0;  // Fermi energy
    double edge   = 0.0;  // x-ray frequency at Fermi level
    double rs     = 0.0;
    double vint   = 0.0;

    // Energy arrays
    std::vector<FeffComplex> em;     // (ne)
    std::vector<FeffComplex> eref;   // (ne, nsp)

    // Per-potential info
    std::array<int, nphx + 1>           iz{};
    std::array<std::string, nphx + 1>   potlbl{};  // 6-char labels

    // Phase shifts  ph(ne, -ltot:ltot, nsp, 0:nph)
    std::vector<FeffComplex> ph;

    // Multipole matrix elements  rkk(ne, 8, nsp)
    std::vector<FeffComplex> rkk;

    // Effective lmax per energy and potential  lmax(ne, 0:nph)
    std::vector<int> lmax;

    // Largest lmax+1 across all energies/potentials
    int lmaxp1 = 0;
};

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Read pot.pad (replaces Fortran rdpot).
PotData read_pot(const std::string& filename = "pot.pad");

/// Read phase.pad (replaces Fortran rdxsph).
PhaseData read_xsph(const std::string& filename = "phase.pad");

/// Read header lines terminated by a "--------" marker (replaces rdhead).
/// Reads from an already-open stream.  Returns lines in `headers`.
void read_header(std::istream& in, std::vector<std::string>& headers);

/// Write header lines (replaces wthead).
void write_header(std::ostream& out, const std::vector<std::string>& headers);

/// Skip comment lines whose first non-whitespace character is in
/// comment_chars (replaces rdcmt).  Leaves the stream positioned at
/// the first non-comment line.
void skip_comments(std::istream& in, const std::string& comment_chars = ";*%#");

/// Build the standard FEFF header block (replaces sthead).
void make_header(int& ntitle, std::array<std::string, nheadx>& title,
                 int nph,
                 const int* iz, const double* rmt, const double* rnrm,
                 const double* xion, int ihole, int ixc,
                 double vr0, double vi0, double gamach,
                 double xmu, double xf, double vint, double rs,
                 int lreal, double rgrd);

} // namespace feff::common
