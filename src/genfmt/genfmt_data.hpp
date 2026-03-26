#pragma once

// GENFMT data structures — C++ equivalents of the Fortran COMMON blocks:
//   /pdata/, /clmz/, /fmatrx/, /lambda/, /nlm/, /rotmat/
// Converted from GENFMT/*.h include files.

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>
#include <array>
#include <complex>
#include <string>

namespace feff::genfmt {

// ---- Dimension constants derived from dim.h ----
// Fortran arrays are 1-based; C++ arrays are 0-based.
// We allocate the same total sizes but shift indexing.

// ---- PhaseData — replaces COMMON /pdata/ ----
// Phase shift data read from phase.pad, shared across genfmt routines.
struct PhaseData {
    // Complex phase shifts: ph(nex, -ltot:ltot, 0:nphx)
    // Stored as [ie][il + ltot][iph] (il offset by ltot for non-negative index)
    FeffComplex ph[nex][2 * ltot + 1][nphx + 1]{};

    // Complex energy reference
    FeffComplex eref[nex]{};

    // Energy mesh
    FeffComplex em[nex]{};

    // Atom positions: rat(3, 0:legtot+1)
    double rat[3][legtot + 2]{};

    // Leg lengths, beta angles, eta angles
    double ri[legtot]{};
    double beta[legtot + 1]{};
    double eta[legtot + 2]{};  // eta(0:legtot+1)

    // Scalar quantities
    double deg = 0.0;
    double rnrmav = 0.0;
    double xmu = 0.0;
    double edge = 0.0;

    // Max l for each potential at each energy: lmax(nex, 0:nphx)
    int lmax[nex][nphx + 1]{};

    // Potential for each atom in path: ipot(0:legtot)
    int ipot[legtot + 1]{};

    // Atomic numbers: iz(0:nphx)
    int iz[nphx + 1]{};

    // Text labels
    std::string text[5];
    int ltext[5]{};
    std::string potlbl[nphx + 1];

    // Scalar integers
    int nsc = 0;
    int nleg = 0;
    int npot = 0;
    int ne = 0;
    int ik0 = 0;
    int ipath = 0;
    int ihole = 0;
    int kinit = 0;
    int linit = 0;
    int ilinit = 0;
    int lmaxp1 = 0;
    int ntext = 0;
};

// ---- LambdaData — replaces COMMON /lambda/ ----
struct LambdaData {
    int mlam[lamtot]{};    // mu for each lambda
    int nlam[lamtot]{};    // nu for each lambda
    int lamx = 0;          // max lambda in problem
    int laml0x = 0;        // max lambda for vectors involving absorbing atom
    int mmaxp1 = 0;        // max mu in problem + 1
    int nmax = 0;          // max nu in problem
};

// ---- NlmData — replaces COMMON /nlm/ ----
struct NlmData {
    // xnlm(ltot+1, mtot+1)
    double xnlm[ltot + 1][mtot + 1]{};
};

// ---- RotationMatrixData — replaces COMMON /rotmat/ ----
struct RotationMatrixData {
    // dri(ltot+1, 2*mtot+1, 2*mtot+1, legtot+1)
    double dri[ltot + 1][2 * mtot + 1][2 * mtot + 1][legtot + 1]{};
};

// ---- ClmzData — replaces COMMON /clmz/ ----
struct ClmzData {
    // clmi(ltot+1, mtot+ntot+1, legtot)
    FeffComplex clmi[ltot + 1][mtot + ntot + 1][legtot]{};
};

// ---- FmatrixData — replaces COMMON /fmatrx/ ----
struct FmatrixData {
    // fmati(lamtot, lamtot, legtot)
    FeffComplex fmati[lamtot][lamtot][legtot]{};
};

// ---- Spin-resolved phase data (4D arrays used in genfmt/onepath) ----
struct SpinPhaseData {
    // ph4(nex, -ltot:ltot, nspx, 0:nphx)
    FeffComplex ph4[nex][2 * ltot + 1][nspx][nphx + 1]{};

    // eref2(nex, nspx)
    FeffComplex eref2[nex][nspx]{};

    // rkk2(nex, 8, nspx)
    FeffComplex rkk2[nex][8][nspx]{};
};

// ---- BmatrixData — used by mmtr/mmtrxi ----
struct BmatrixData {
    // bmati(-mtot:mtot, 8, -mtot:mtot, 8)
    // stored as [m1 + mtot][k1][m2 + mtot][k2]
    FeffComplex bmati[2 * mtot + 1][8][2 * mtot + 1][8]{};
};

// Indexing helpers for Fortran-style shifted arrays.
// Phase shift: ph(ie, il, iph) where il in [-ltot, ltot]
inline FeffComplex& ph_at(PhaseData& pd, int ie, int il, int iph) {
    return pd.ph[ie][il + ltot][iph];
}
inline const FeffComplex& ph_at(const PhaseData& pd, int ie, int il, int iph) {
    return pd.ph[ie][il + ltot][iph];
}

// Spin phase: ph4(ie, il, is, iph) where il in [-ltot, ltot], is 0-based
inline FeffComplex& ph4_at(SpinPhaseData& sp, int ie, int il, int is, int iph) {
    return sp.ph4[ie][il + ltot][is][iph];
}
inline const FeffComplex& ph4_at(const SpinPhaseData& sp, int ie, int il, int is, int iph) {
    return sp.ph4[ie][il + ltot][is][iph];
}

// Rotation matrix: dri(il, m1, m2, ileg) where m1,m2 in [-mtot, mtot]
inline double& dri_at(RotationMatrixData& rm, int il, int m1, int m2, int ileg) {
    return rm.dri[il][m1 + mtot][m2 + mtot][ileg];
}
inline double dri_at(const RotationMatrixData& rm, int il, int m1, int m2, int ileg) {
    return rm.dri[il][m1 + mtot][m2 + mtot][ileg];
}

// Bmati: bmati(m1, k1, m2, k2) where m1,m2 in [-mtot, mtot], k1,k2 in [0,7]
inline FeffComplex& bmati_at(BmatrixData& bd, int m1, int k1, int m2, int k2) {
    return bd.bmati[m1 + mtot][k1][m2 + mtot][k2];
}
inline const FeffComplex& bmati_at(const BmatrixData& bd, int m1, int k1, int m2, int k2) {
    return bd.bmati[m1 + mtot][k1][m2 + mtot][k2];
}

// eta: eta(j) where j in [0, legtot+1]
inline double& eta_at(PhaseData& pd, int j) {
    return pd.eta[j];
}
inline double eta_at(const PhaseData& pd, int j) {
    return pd.eta[j];
}

} // namespace feff::genfmt
