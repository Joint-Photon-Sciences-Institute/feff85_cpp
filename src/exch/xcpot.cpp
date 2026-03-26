// Main exchange-correlation potential driver
// Converted from src/EXCH/xcpot.f
//
// First coded by J. Mustre de Leon
// Last modified by A. Ankudinov 1996 for non-local self-energies
// Many-pole additions by Josh Kas

#include "xcpot.hpp"
#include "csigma.hpp"
#include "csigz.hpp"
#include "edp.hpp"
#include "vbh.hpp"
#include "rhl.hpp"
#include "rhlbp.hpp"
#include "imhl.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include <feff/types.hpp>
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../math/interpolation.hpp"
#include "../common/logging.hpp"
#include "../par/parallel.hpp"

namespace feff::exch {

// ---------------------------------------------------------------------------
// Internal sigma dispatcher (from the nested subroutine sigma in xcpot.f)
// ---------------------------------------------------------------------------
static void sigma_dispatch(int ixc, int ibp, double rs, double rscore,
                           double xk, double& vr, double& vi)
{
    vr = 0.0;
    vi = 0.0;

    if ((ixc == 0 || ixc >= 5) && ibp == 0) {
        rhl(rs, xk, vr, vi);
    } else if ((ixc == 0 || ixc >= 5) && ibp == 1) {
        rhlbp(rs, xk, vr, vi);
    } else if (ixc == 1) {
        vi = 0.0;
        edp(rs, xk, vr);
    } else if (ixc == 3) {
        edp(rs, xk, vr);
        int icusp = 0;
        imhl(rs, xk, vi, icusp);
    }

    if (ixc >= 6) {
        double vrp = 0.0;
        edp(rscore, xk, vrp);
        vr = vr - vrp;
    }
}

// ---------------------------------------------------------------------------
// xcpot main routine
// ---------------------------------------------------------------------------
void xcpot(int iph, int ie, int index, int lreal, int& ifirst, int jri,
           FeffComplex em, double xmu,
           const double vtot[], const double vvalgs[],
           const double densty[], const double dmag[], const double denval[],
           FeffComplex& eref, FeffComplex v[], FeffComplex vval[],
           int ipl, const double wpcorr[], const double ampfac[],
           double vxcrmu[], double vxcimu[], double gsrel[],
           double vvxcrm[], double vvxcim[], double rnrm)
{
    constexpr int NRPts = 10;
    constexpr double dx = 0.05;
    constexpr double x0 = 8.8;

    double totvol = 4.0 / 3.0 * pi * rnrm * rnrm * rnrm;
    FeffComplex delavg(0.0, 0.0);
    FeffComplex zrnrm(0.0, 0.0);

    bool csig = false;
    int lastpl = 0;
    int nmax = 1;
    int nul = 0;

    int ibp = index / 10;
    int ixc = index % 10;

    // Find last valid pole
    for (int i = 0; i < MxPole; ++i) {
        if (wpcorr[i] <= 0.0) {
            lastpl = i;  // 0-based: lastpl is count of valid poles
            break;
        }
    }

    if (ixc == 0 && ipl > 0) {
        csig = true;
    }

    // Ground state exchange or below Fermi level: no self-energy
    if (ixc == 2 || em.real() <= xmu) {
        for (int i = 0; i <= jri; ++i) {
            v[i] = vtot[i];
            vval[i] = vvalgs[i];
        }
        // Skip to reference potential (goto 888)
        eref = v[jri];
        for (int i = 0; i <= jri; ++i) {
            v[i] = v[i] - eref;
        }
        if (ixc >= 5) {
            for (int i = 0; i <= jri; ++i) {
                vval[i] = vval[i] - eref;
            }
        } else {
            for (int i = 0; i <= jri; ++i) {
                vval[i] = v[i];
            }
        }
        if (lreal > 0) {
            for (int i = 0; i <= jri; ++i) {
                v[i] = v[i].real();
                if (ixc > 4) vval[i] = vval[i].real();
            }
            eref = eref.real();
        }
        // NOTE: Do NOT set ifirst=1 here. In the Fortran original, the
        // "goto 888" jumps past "ifirst = 1", so when the ground-state /
        // below-Fermi early return is taken, ifirst must remain 0.
        // This ensures that on the next call (with em > xmu), the
        // Fermi-level self-energy arrays vxcrmu/vxcimu are computed.
        return;
    }

    // Calculate Rs at core and interstitial densities
    double rsint;
    if (densty[jri] <= 0.0) {
        rsint = 10.0;
    } else {
        rsint = std::pow(3.0 / (4.0 * pi * densty[jri]), third);
    }

    double rscore;
    if (densty[0] <= 0.0) {
        rscore = 101.0;
    } else {
        rscore = std::pow(3.0 / (4.0 * pi * densty[0]), third);
    }

    double drs = (rsint - rscore) / (NRPts - 1);
    double omp_int = std::sqrt(3.0 / (rsint * rsint * rsint)) * hart;
    double ompmax = omp_int * wpcorr[lastpl > 0 ? lastpl - 1 : 0];

    // Calculate delta sigma as a function of Rs and energy for many-pole model
    double delrHL[NRPts] = {};
    double deliHL[NRPts] = {};
    double rs1[NRPts] = {};

    if (csig) {
        for (int i = NRPts - 1; i >= 0; --i) {
            delrHL[i] = 0.0;
            deliHL[i] = 0.0;
            rs1[i] = rscore + static_cast<double>(i) * drs;

            FeffComplex ztemp(0.0, 0.0);

            if (ipl > 1) {
                if (ipl == 2 || i == NRPts - 1) {
                    csigz(em, xmu, rs1[i], delrHL[i], deliHL[i], ztemp,
                          wpcorr, ampfac);
                } else if (ipl == 3) {
                    delrHL[i] = 0.0;
                    deliHL[i] = 0.0;
                } else {
                    delrHL[i] = delrHL[NRPts - 1];
                    deliHL[i] = deliHL[NRPts - 1];  // Fortran: only delrHL was assigned
                }
                if (i == NRPts - 1) zrnrm = ztemp;
            } else {
                csigma(em, xmu, rs1[i], delrHL[i], deliHL[i], wpcorr, ampfac);
            }
        }
    }

    // Add the self energy correction, loop from jri1 down to 0
    // (Fortran loop: do 20 i = jri1, 1, -1 => i is 1-based)
    double riint = 0.0;
    double xfval = 0.0;
    double delvr = 0.0;
    double delvi = 0.0;

    for (int i = jri; i >= 0; --i) {
        double ri = std::exp(static_cast<double>(i) * dx - x0);
        if (i == jri) {
            riint = ri;
        }
        int niter = 0;

        double rs_local;
        if (densty[i] <= 0.0) {
            rs_local = 10.0;
        } else {
            rs_local = std::pow(3.0 / (4.0 * pi * densty[i]), third);
        }

        double delr = 0.0, deli = 0.0;

        // If csigma is on, interpolate onto rs
        if (csig) {
            double omp_local = std::sqrt(3.0 / (rs_local * rs_local * rs_local)) * hart;
            if (ipl >= 4) {
                delr = delrHL[NRPts - 1];
                deli = deliHL[NRPts - 1];
            } else {
                feff::math::terp(rs1, delrHL, NRPts, 1, rs_local, delr);
                feff::math::terp(rs1, deliHL, NRPts, 1, rs_local, deli);
            }
            if (ipl != 5 || omp_local < ompmax) {
                // Skip SC loop, go directly to forming delta
                FeffComplex delta(delr, deli);
                if (ixc == 5) delta = FeffComplex(delr, delvi);
                v[i] = vtot[i] + delta;
                if (ixc >= 5) {
                    FeffComplex deltav(delvr, delvi);
                    vval[i] = vvalgs[i] + deltav;
                }

                // Volume calculation for averaging
                double volume;
                if (i == jri) {
                    volume = 0.0;
                } else if (i == jri - 1) {
                    volume = 4.0 / 3.0 * pi * (rnrm * rnrm * rnrm
                             - std::exp(3.0 * (static_cast<double>(i) * dx - x0)));
                } else {
                    volume = 4.0 / 3.0 * pi * std::exp(3.0 * (static_cast<double>(i) * dx - x0))
                             * (std::exp(3.0 * dx) - 1.0);
                }
                if (volume < 0.0) volume = 0.0;

                continue;
            }
        }

        // Standard self-energy calculation
        double xf = fa / rs_local;
        double rsm = rs_local / std::pow(1.0 + dmag[i], third);
        double xfm = fa / rsm;

        double rsval = 0.0;
        if (ixc == 5) {
            if (denval[i] > 0.00001) {
                rsval = std::pow(3.0 / (4.0 * pi * denval[i]), third);
                if (rsval > 10.0) rsval = 10.0;
            } else {
                rsval = 10.0;
            }
            xfval = fa / rsval;
        } else if (ixc >= 6) {
            if (densty[i] <= denval[i]) {
                rscore = 101.0;
            } else {
                rscore = std::pow(3.0 / (4.0 * pi * (densty[i] - denval[i])), third);
            }
        }

        if (ifirst == 0) {
            // vxc at Fermi level, calculate only once
            double xk = xf * 1.00001;
            gsrel[i] = 1.0;
            if (ixc < 5) {
                sigma_dispatch(ixc, ibp, rs_local, rscore, xk, vxcrmu[i], vxcimu[i]);
            } else {
                sigma_dispatch(nul, ibp, rs_local, rscore, xk, vxcrmu[i], vxcimu[i]);
            }

            if (ixc == 5) {
                double xkpp = xfval * 1.00001;
                sigma_dispatch(ixc, ibp, rsval, rscore, xkpp, vvxcrm[i], vvxcim[i]);
                if (ixc == 5 && i == jri) {
                    vvxcrm[jri] = vxcrmu[jri];
                    vvxcim[jri] = vxcimu[jri];
                }
            } else if (ixc >= 6) {
                sigma_dispatch(ixc, ibp, rs_local, rscore, xk, vvxcrm[i], vvxcim[i]);
                if (ixc == 6 && i == jri) {
                    vvxcrm[jri] = vxcrmu[jri];
                    vvxcim[jri] = vxcimu[jri];
                }
            } else {
                vvxcrm[i] = 0.0;
                vvxcim[i] = 0.0;
            }
        }

        // xk2 is the local momentum squared: p^2 = k^2 - 2*mu + kf^2
        double xk2 = 2.0 * (em.real() - xmu) + xf * xf;
        double xk = std::sqrt(xk2);
        double xkm2 = 2.0 * (em.real() - xmu) + xfm * xfm;
        if (xkm2 < 0.0) xkm2 = xk2;
        double xkm = std::sqrt(xkm2);

        // Find delta_1
        double vxcr = 0.0, vxci = 0.0;
        if (ixc < 5) {
            sigma_dispatch(ixc, ibp, rs_local, rscore, xk, vxcr, vxci);
        } else {
            sigma_dispatch(nul, ibp, rs_local, rscore, xk, vxcr, vxci);
        }
        double del1r = gsrel[i] * (vxcr - vxcrmu[i]);

        // Iterative solution of Dyson equation
        for (;;) {
            xk2 = 2.0 * (em.real() - xmu - del1r) + xf * xf;
            if (xk2 < 0.0) {
                std::ostringstream msg;
                msg << xk2 << " " << i << " " << ie << " " << iph
                    << " xk2, i, ie, iph";
                feff::common::logger().wlog(msg.str());
                feff::common::logger().wlog(" em, xf**2, xmu, delta");
                msg.str("");
                msg << em.real() << " " << xf * xf << " " << xmu << " " << del1r;
                feff::common::logger().wlog(msg.str());
                feff::par::par_stop("XCPOT-2");
            }
            xk = std::sqrt(xk2);

            // Calculate delta_2 with corrected local momentum
            sigma_dispatch(ixc, ibp, rs_local, rscore, xk, vxcr, vxci);
            delr = gsrel[i] * (vxcr - vxcrmu[i]);
            deli = vxci - vxcimu[i];

            if (ixc >= 5 && i == jri && xk > xf) {
                if (ixc == 5 || ixc == 6) {
                    delvr = delr;
                    delvi = deli;
                }
            }

            if (niter < nmax) {
                del1r = delr;
                ++niter;
                continue;
            }
            break;
        }

        if (ixc >= 5 && i < jri && xk > xf) {
            double vxcvr = 0.0, vxcvi = 0.0;
            if (ixc == 5) {
                double xkpp = std::sqrt(xk * xk - xf * xf + xfval * xfval);
                sigma_dispatch(ixc, ibp, rsval, rscore, xkpp, vxcvr, vxcvi);
            } else {
                sigma_dispatch(ixc, ibp, rs_local, rscore, xk, vxcvr, vxcvi);
            }
            delvr = vxcvr - vvxcrm[i];
            delvi = vxcvi - vvxcim[i];
        }

        // Form complex delta and potential
        FeffComplex delta(delr, deli);
        if (ixc == 5) delta = FeffComplex(delr, delvi);

        v[i] = vtot[i] + delta;

        if (ixc >= 5) {
            FeffComplex deltav(delvr, delvi);
            vval[i] = vvalgs[i] + deltav;
        }

        // Volume calculation
        double volume;
        if (i == jri) {
            volume = 0.0;
        } else if (i == jri - 1) {
            volume = 4.0 / 3.0 * pi * (rnrm * rnrm * rnrm
                     - std::exp(3.0 * (static_cast<double>(i) * dx - x0)));
        } else {
            volume = 4.0 / 3.0 * pi * std::exp(3.0 * (static_cast<double>(i) * dx - x0))
                     * (std::exp(3.0 * dx) - 1.0);
        }
        if (volume < 0.0) volume = 0.0;
    }

    ifirst = 1;
    delavg = delavg / totvol;

    // Reference potential with respect to MT potential (v[jri] = 0)
    eref = v[jri];
    for (int i = 0; i <= jri; ++i) {
        v[i] = v[i] - eref;
    }
    if (ixc >= 5) {
        for (int i = 0; i <= jri; ++i) {
            vval[i] = vval[i] - eref;
        }
    } else {
        for (int i = 0; i <= jri; ++i) {
            vval[i] = v[i];
        }
    }

    // Real self energy: zero imag part
    if (lreal > 0) {
        for (int i = 0; i <= jri; ++i) {
            v[i] = v[i].real();
            if (ixc > 4) vval[i] = vval[i].real();
        }
        eref = eref.real();
    }
}

} // namespace feff::exch
