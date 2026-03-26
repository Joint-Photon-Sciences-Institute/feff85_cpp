#pragma once

// FeffInput: aggregate struct holding all parsed FEFF input data
// Converted from Fortran common blocks in RDINP/allinp.h

#include <array>
#include <complex>
#include <string>
#include <vector>

#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff {

// Version strings (from HEADERS/vers.h)
inline const std::string vfeff_default = "Feff8L (EXAFS)";
inline const std::string vf85e_default = " 0.1";

struct FeffInput {
    // --- Version info ---
    std::string vfeff = vfeff_default;
    std::string vf85e = vf85e_default;

    // --- ATOMS (atoms.dat) --- common /geom/
    int natt = 0;
    std::array<int, nattx> iphatx{};          // iphatx(nattx)
    std::array<Vec3, nattx> ratx{};           // ratx(3, nattx) stored as array of Vec3

    // --- GLOBAL (global.inp) --- common /global/
    int iphabs = 0;
    int ipol = 0;
    int ispin = 0;
    int le2 = 0;
    Vec3 evec{};
    Vec3 xivec{};
    Vec3 spvec{};
    double elpty = 0.0;
    double angks = 0.0;
    double rclabs = 0.0;
    // ptz(-1:1, -1:1) -> 3x3 complex array, indexed [row+1][col+1]
    std::array<std::array<FeffComplex, 3>, 3> ptz{};

    // --- MOD1 (pot.json) --- common /mod1/
    std::array<std::string, nheadx> title{};
    int mpot = 0;
    int nph = 0;
    int ntitle = 0;
    int ihole = 0;
    int ipr1 = 0;
    int iafolp = 0;
    int iunf = 0;
    int nmix = 0;
    int nohole = 0;
    int jumprm = 0;
    int inters = 0;
    int nscmt = 0;
    int icoul = 0;
    int lfms1 = 0;

    std::array<int, nphx + 1> iz{};           // iz(0:nphx)
    std::array<int, nphx + 1> lmaxsc{};       // lmaxsc(0:nphx)
    float rfms1 = 0.0f;
    double gamach = 0.0;
    double rgrd = 0.0;
    double ca1 = 0.0;
    double ecv = 0.0;
    double totvol = 0.0;
    std::array<double, nphx + 1> xnatph{};    // xnatph(0:nphx)
    std::array<double, nphx + 1> folp{};      // folp(0:nphx)
    std::array<double, nphx + 1> spinph{};    // spinph(0:nphx)
    std::array<double, nphx + 1> xion{};      // xion(0:nphx)

    // OVERLAP arrays
    std::array<int, nphx + 1> novr{};                         // novr(0:nphx)
    std::array<std::array<int, novrx>, nphx + 1> iphovr{};    // iphovr(novrx, 0:nphx)
    std::array<std::array<int, novrx>, nphx + 1> nnovr{};     // nnovr(novrx, 0:nphx)
    std::array<std::array<double, novrx>, nphx + 1> rovr{};   // rovr(novrx, 0:nphx)

    // --- MOD2 (xsph.json) --- common /mod2/
    int mphase = 0;
    int ipr2 = 0;
    int ixc = 0;
    int ixc0 = 0;
    int ispec = 0;
    int lreal = 0;
    int l2lp = 0;
    int iPlsmn = 0;
    int iGrid = 0;
    int lfms2 = 0;
    double vr0 = 0.0;
    double vi0 = 0.0;
    double rfms2 = 0.0;
    double xkstep = 0.0;
    double xkmax = 0.0;
    double vixan = 0.0;
    int izstd = 0;
    int ifxc = 0;
    int ipmbse = 0;
    int itdlda = 0;
    int nonlocal = 0;
    int ibasis = 0;
    std::array<int, nphx + 1> lmaxph{};       // lmaxph(0:nphx)
    std::array<std::string, nphx + 1> potlbl{};  // potlbl(0:nphx), char*6

    // --- MOD3 (fms.json) --- common /mod3/
    int mfms = 0;
    int idwopt = 0;
    int minv = 0;
    float rdirec = 0.0f;
    float toler1 = 0.0f;
    float toler2 = 0.0f;
    double tk = 0.0;
    double thetad = 0.0;
    double sig2g = 0.0;

    // --- MOD4 (path.json) --- common /mod4/
    int mpath = 0;
    int ms = 0;
    int nncrit = 0;
    int nlegxx = 0;
    int ipr4 = 0;
    float critpw = 0.0f;
    float pcritk = 0.0f;
    float pcrith = 0.0f;
    float rmax = 0.0f;

    // --- MOD5 (genfmt.json) --- common /mod5/
    int mfeff = 0;
    int ipr5 = 0;
    int iorder = 0;
    bool wnstar = false;
    double critcw = 0.0;

    // --- MOD6 (ff2x.json) --- common /mod6/
    int mchi = 0;
    int ipr6 = 0;
    int mbconv = 0;
    int absolu = 0;
    double vrcorr = 0.0;
    double vicorr = 0.0;
    double s02 = 0.0;
    double alphat = 0.0;
    double thetae = 0.0;

    // --- SO2 (s02.json) --- common /so2/
    int mso2conv = 0;
    int ipse = 0;
    int ipsk = 0;
    double wsigk = 0.0;
    double cen = 0.0;
    std::string cfname;

    // --- LDOS / energy grid ---
    int mldos = 0;
    double emin = 0.0;
    double emax = 0.0;
    double eimag = 0.0;

    // --- EELS (common /eelsva/) ---
    int eels = 0;
    int relat = 0;
    int aver = 0;
    int cross = 0;
    int iinput = 0;
    int spcol = 0;
    int nqr = 0;
    int nqf = 0;
    int magic = 0;
    int ipmin = 0;
    int ipmax = 0;
    int ipstep = 0;
    double ebeam = 0.0;
    double aconv = 0.0;
    double acoll = 0.0;
    double thetax = 0.0;
    double thetay = 0.0;
    double emagic = 0.0;
};

} // namespace feff
