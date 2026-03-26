#pragma once
// Master orchestrator for phase shifts and cross-sections.
// Converted from src/XSPH/xsph.f (~570 lines)
//
// This is the main entry point for the XSPH module. It:
// 1. Constructs the energy mesh
// 2. Calculates the cross-section for the absorbing atom
// 3. Calculates phase shifts for all unique potentials
// 4. Writes output files (phase.pad, xsect.json, axafs.dat)

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include <string>

namespace feff::xsph {

/// Main XSPH calculation: phase shifts and cross-sections.
///
/// @param wrxsec  Write xsect.json output
/// @param verbose Enable log messages
/// @param phpad   Output filename for phase.pad
/// @param ipr2    Print level
/// @param ispec   Calculation type (0=EXAFS, 1-3=XANES/XES/DANES, 4=FPRIME)
/// @param vixan   FMS energy step (Hartrees)
/// @param xkstep  k-step (inverse Bohr)
/// @param xkmax   Maximum k (inverse Bohr)
/// @param gamach  Core-hole broadening (Hartrees)
/// @param rgrd    Phase grid step
/// @param nph     Number of unique potentials
/// @param lmaxph  Max l per potential [nphx+1]
/// @param potlbl  Potential labels [nphx+1][7]
/// @param spinph  Spin density factors [nphx+1]
/// @param iatph   Representative atom per potential [nphx+1]
/// @param nat     Number of atoms
/// @param rat     Atomic coordinates [3][natx]
/// @param iphat   Potential assignment per atom [natx]
/// @param ixc     Exchange-correlation model
/// @param vr0     Real part of constant potential shift (Hartrees)
/// @param vi0     Imaginary part of constant potential (Hartrees)
/// @param ixc0    XC model for cross-section
/// @param lreal   Real phase shift flag
/// @param rfms2   FMS radius
/// @param lfms2   FMS angular momentum flag
/// @param l2lp    Transition selector
/// @param ipol    Polarization type
/// @param ispin   Spin flag
/// @param le2     Multipole selector
/// @param angks   k-spin angle
/// @param ptz     Polarization tensor [3][3]
/// @param iPl     Many-pole self-energy control
/// @param izstd   TDLDA flag
/// @param ifxc    Core-hole interaction flag
/// @param ipmbse  Many-body flag
/// @param itdlda  TDLDA iteration flag
/// @param nonlocal Nonlocal potential flag
/// @param ntitle  Number of title lines
/// @param title   Title strings [nheadx][80]
/// @param rnrmav  Average Norman radius
/// @param xmu     Fermi level (Hartrees)
/// @param vint    Interstitial potential (Hartrees)
/// @param rhoint  Interstitial density
/// @param emu     Edge energy (Hartrees)
/// @param s02     S0^2 amplitude factor
/// @param erelax  Relaxation energy
/// @param wp      Plasmon frequency
/// @param ecv     Core-valence separation
/// @param rs      Density parameter
/// @param xf      Fermi momentum
/// @param qtotel  Total charge
/// @param imt     MT grid indices [nphx+1]
/// @param rmt     MT radii [nphx+1]
/// @param inrm    Norman grid indices [nphx+1]
/// @param rnrm    Norman radii [nphx+1]
/// @param folp    Overlap fractions [nphx+1]
/// @param folpx   Max overlap fractions [nphx+1]
/// @param xnatph  Number of atoms per potential [nphx+1]
/// @param dgc0    Initial orbital large component [251]
/// @param dpc0    Initial orbital small component [251]
/// @param dgc     All large Dirac components [251][30][nphx+2]
/// @param dpc     All small Dirac components [251][30][nphx+2]
/// @param adgc    Development coefficients [10][30][nphx+2]
/// @param adpc    Development coefficients [10][30][nphx+2]
/// @param edens   Electron density [251][nphx+1]
/// @param vclap   Coulomb potential [251][nphx+1]
/// @param vtot    Total potential [251][nphx+1]
/// @param edenvl  Valence density [251][nphx+1]
/// @param vvalgs  Valence potential [251][nphx+1]
/// @param dmag    Magnetization density [251][nphx+2]
/// @param xnval   Valence occupations [30][nphx+1]
/// @param iorb    Orbital indices [-4:3][nphx+1]
/// @param nohole  No-hole flag
/// @param ihole   Core-hole orbital index
/// @param inters  Interstitial model flag
/// @param totvol  Total volume
/// @param iafolp  Automatic FOLP flag
/// @param xion    Ionicities [nphx+1]
/// @param iunf    Unfolding flag
/// @param iz      Atomic numbers [nphx+1]
/// @param jumprm  Jump flag for potential
void xsph(bool wrxsec, bool verbose, const std::string& phpad,
           int ipr2, int ispec, double vixan, double xkstep, double xkmax,
           double gamach, double rgrd,
           int nph, int lmaxph[], const char potlbl[][7],
           const double spinph[], const int iatph[], int nat,
           const double* rat, const int iphat[],
           int ixc, double vr0, double vi0, int ixc0, int lreal,
           float rfms2, int lfms2, int l2lp,
           int ipol, int ispin, int le2, double angks, const FeffComplex ptz[3][3],
           int iPl, int izstd, int ifxc, int ipmbse, int itdlda, int nonlocal,
           int ntitle, const char title[][80], double rnrmav,
           double xmu, double vint, double rhoint,
           double emu, double s02, double erelax, double wp, double ecv,
           double rs, double xf, double qtotel,
           const int imt[], const double rmt[], const int inrm[],
           const double rnrm[], const double folp[], double folpx[],
           const double xnatph[],
           const double dgc0[], const double dpc0[],
           const double* dgc, const double* dpc,
           const double* adgc, const double* adpc,
           double* edens, double* vclap, double* vtot,
           double* edenvl, double* vvalgs, double* dmag,
           const double* xnval, const int* iorb,
           int nohole, int ihole,
           int inters, double totvol, int iafolp,
           const double xion[], int iunf, const int iz[], int jumprm);

} // namespace feff::xsph
