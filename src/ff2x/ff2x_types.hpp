#pragma once
// Data types and structures for the FF2X module.
// Converted from Fortran FF2X module types and common blocks.

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>

#include <array>
#include <complex>
#include <string>
#include <vector>

namespace feff::ff2x {

// Maximum number of paths in list.dat
inline constexpr int npx = 15000;

// Maximum number of fine-grid points for EXAFS chi(k)
inline constexpr int nfinex = 601;

// Small numerical tolerance
inline constexpr double eps4 = 1.0e-4;

// Comment prefix for output files
inline constexpr const char* coment = "# ";

// ---------------------------------------------------------------------------
// FeffPadData -- data read from feff.pad (PAD-format binary)
// ---------------------------------------------------------------------------
struct FeffPadData {
    int nptot  = 0;   // total number of paths
    int ne     = 0;   // number of energy points
    int npot   = 0;   // number of potential types
    int ihole  = 0;
    int iorder = 0;
    int ilinit = 0;

    float rnrmav = 0.0f;   // average Norman radius
    float xmu    = 0.0f;   // Fermi level
    float edge   = 0.0f;   // edge energy (code units)

    std::array<std::string, nphx + 1> potlbl{};  // potential labels (0:nphx)
    std::array<int, nphx + 1> iz{};               // atomic numbers (0:nphx)

    // Central atom phase shift (complex single precision stored as double)
    std::array<std::complex<float>, nex> phc{};
    std::array<std::complex<float>, nex> ck{};
    std::array<float, nex> xk{};

    // Per-path arrays (up to npx paths)
    std::vector<int> index;               // (npx) path index
    std::vector<int> nleg;                // (npx) number of legs
    std::vector<float> deg;               // (npx) degeneracy
    std::vector<float> reff;              // (npx) half path length
    std::vector<float> crit;              // (npx) importance factor
    std::vector<std::vector<int>> ipot;   // (legtot, npx) potential indices
    std::vector<std::vector<std::array<float, 3>>> rat;  // (3, legtot, npx) positions
    std::vector<std::vector<float>> beta; // (legtot, npx)
    std::vector<std::vector<float>> eta;  // (legtot, npx)
    std::vector<std::vector<float>> ri;   // (legtot, npx)

    // Per-path amplitude & phase: achi(nex, npx), phchi(nex, npx)
    std::vector<std::vector<float>> achi;
    std::vector<std::vector<float>> phchi;
};

// ---------------------------------------------------------------------------
// XsectData -- data read from xsect.bin (via rdxbin / read_xsect)
// ---------------------------------------------------------------------------
struct XsectData {
    double s02p    = 0.0;
    double s02_eff = 1.0;  // effective s02: overridden from s02p when input s02 <= 0.1
    double erelax  = 0.0;
    double wp     = 0.0;
    double edgep  = 0.0;
    double gamach  = 0.0;
    int ne1       = 0;
    int ik0       = 0;
    int nxsec     = 0;
    int mbconv    = 0;
    int ntitle    = 0;

    std::array<std::string, nheadx> title{};
    std::array<FeffComplex, nex> emxs{};   // complex energy grid
    std::array<FeffComplex, nex> xsec{};   // cross section
    std::array<double, nex> omega{};       // energy / work array
    std::array<double, nex> xkxs{};        // k-grid from xsect
    std::array<double, nex> xsnorm{};      // normalized cross-section
};

// ---------------------------------------------------------------------------
// FF2xParams -- input parameters for ff2x (from ff2x.json + global.json)
// ---------------------------------------------------------------------------
struct FF2xParams {
    int mchi    = 0;     // calculation flag
    int ispec   = 0;     // spectroscopy type (0=EXAFS, 1=XANES, 2=XES, 3=DANES, 4=FPRIME)
    int ipr6    = 0;     // print level
    int idwopt  = 0;     // DW option (0=CD, 1=EM, 2=RM, 3=CL)
    int mbconv  = 0;     // convolution flag
    int absolu  = 0;     // absolute normalization flag

    double vrcorr  = 0.0;  // real energy correction (eV, converted to code units)
    double vicorr  = 0.0;  // imaginary energy correction (eV, converted to code units)
    double s02     = 0.0;  // amplitude reduction factor
    double critcw  = 0.0;  // curved-wave amplitude filter (%)
    double tk      = 0.0;  // temperature (K)
    double thetad  = 0.0;  // Debye temperature (K)
    double alphat  = 0.0;  // thermal expansion coefficient
    double thetae  = 0.0;  // Einstein temperature (K)
    double sig2g   = 0.0;  // global sigma^2

    // From global.json
    int nabs    = 1;     // number of absorbers for configuration average
    int iphabs  = 0;
    double rclabs = 0.0;

    // Spectrum loop control
    int ipmin   = 1;
    int ipmax   = 1;
    int ipstep  = 1;
    int elnes   = 0;
};

// ---------------------------------------------------------------------------
// PathListEntry -- one entry from list.dat
// ---------------------------------------------------------------------------
struct PathListEntry {
    int ip = 0;          // path index
    double sig2u = 0.0;  // user-specified sigma^2
};

} // namespace feff::ff2x
