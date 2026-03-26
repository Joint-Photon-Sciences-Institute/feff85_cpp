// Initialize all input variables to their defaults.
// Converted from: src/RDINP/iniall.f

#include "iniall.hpp"
#include <feff/dimensions.hpp>

namespace feff::rdinp {

void iniall(FeffInput& inp) {

    // Initialize title lines
    for (auto& t : inp.title) t = " ";

    // Initialize integer scalars
    inp.iGrid = 0;
    inp.ntitle = 0;
    // nat is not in FeffInput (local to ffsort); natt is set below
    inp.natt = 0;
    inp.nph = 0;

    inp.iafolp = 0;
    inp.idwopt = -1;
    inp.ihole = 1;
    inp.inters = 0;
    inp.iorder = 2;
    inp.ipr1 = 0;
    inp.ipr2 = 0;
    // ipr3 is local to rdinp, not in FeffInput
    inp.ipr4 = 0;
    inp.ipr5 = 0;
    inp.ipr6 = 0;
    inp.ipse = 0;
    inp.ipsk = 0;
    inp.ispec = 0;
    inp.ixc = 0;
    inp.ixc0 = -1;
    inp.jumprm = 0;
    inp.lfms1 = 0;
    inp.lfms2 = 0;
    inp.minv = 0;
    inp.lreal = 0;
    inp.mbconv = 0;
    inp.mchi = 1;
    inp.mfeff = 1;
    inp.mfms = 1;
    inp.mpath = 1;
    inp.mphase = 1;
    inp.mldos = 0;
    inp.mpot = 1;
    inp.ms = 0;
    inp.iPlsmn = 0;
    inp.mso2conv = 0;
    inp.nlegxx = 10;
    inp.nmix = 1;
    inp.nohole = -1;
    inp.nscmt = 0;
    inp.icoul = 0;
    inp.iunf = 0;
    inp.izstd = 0;
    inp.ifxc = 0;
    inp.ipmbse = 0;
    inp.itdlda = 0;
    inp.nonlocal = 0;
    inp.ibasis = 0;

    // Initialize reals (single precision in Fortran)
    inp.critpw = 2.5f;
    inp.pcritk = 0.0f;
    inp.pcrith = 0.0f;
    inp.rmax = -1.0f;
    inp.rfms1 = -1.0f;
    inp.rfms2 = -1.0;
    inp.rdirec = -1.0f;
    inp.toler1 = 1.0e-3f;
    inp.toler2 = 1.0e-3f;

    // Initialize double precision scalars
    inp.alphat = 0.0;
    inp.thetae = 0.0;
    inp.ca1 = 0.0;
    inp.critcw = 4.0;
    inp.eimag = -1.0;
    inp.ecv = -40.0;
    inp.emax = 0.0;
    inp.emin = 1000.0;
    inp.rclabs = 0.0;
    inp.rgrd = 0.05;
    inp.s02 = 1.0;
    inp.sig2g = 0.0;
    inp.thetad = 0.0;
    inp.tk = 0.0;
    inp.totvol = 0.0;
    inp.vr0 = 0.0;
    inp.vi0 = 0.0;
    inp.vicorr = 0.0;
    inp.vrcorr = 0.0;
    inp.xkmax = 20.0;
    inp.xkstep = 0.07;
    inp.vixan = 0.0;
    inp.wsigk = 0.0;
    inp.cen = 0.0;

    // Initialize logicals
    inp.wnstar = false;

    // Initialize per-potential arrays (0:nphx)
    for (int i = 0; i <= nphx; ++i) {
        inp.xnatph[i] = 0.0;
        inp.spinph[i] = 0.0;
        inp.iz[i] = 0;
        inp.xion[i] = 0.0;
        inp.folp[i] = 1.0;
        inp.novr[i] = 0;
        inp.lmaxsc[i] = 0;
        inp.lmaxph[i] = 0;
        inp.potlbl[i] = " ";
    }

    // Initialize overlap arrays (novrx, 0:nphx)
    for (int i = 0; i <= nphx; ++i) {
        for (int j = 0; j < novrx; ++j) {
            inp.iphovr[i][j] = 0;
            inp.nnovr[i][j] = 0;
            inp.rovr[i][j] = 0.0;
        }
    }

    // Initialize polarization data
    inp.ipol = 0;
    inp.ispin = 0;
    inp.le2 = 0;
    inp.l2lp = 0;
    inp.elpty = 0.0;
    inp.angks = 0.0;
    for (int i = 0; i < 3; ++i) {
        inp.evec[i] = 0.0;
        inp.xivec[i] = 0.0;
        inp.spvec[i] = 0.0;
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            inp.ptz[i][j] = FeffComplex(0.0, 0.0);
        }
    }

    // Initialize atom list data
    inp.iphatx.fill(-1);
    inp.ratx.fill(Vec3{0.0, 0.0, 0.0});

    // Initialize character strings
    inp.cfname = "NULL";

    // Initialize EELS variables
    inp.ebeam = 0.0;
    inp.aconv = 0.0;
    inp.acoll = 0.0;
    inp.nqr = 0;
    inp.nqf = 0;
    inp.magic = 0;
    inp.emagic = 0.0;
    inp.eels = 0;
    inp.relat = 1;
    inp.cross = 1;
    inp.aver = 0;
    inp.thetax = 0.0;
    inp.thetay = 0.0;
    inp.ipmin = 1;
    inp.ipmax = 9;
    inp.ipstep = 1;
    inp.iinput = 1;
    inp.spcol = 4;

    // ABSOLUTE card
    inp.absolu = 0;
}

} // namespace feff::rdinp
