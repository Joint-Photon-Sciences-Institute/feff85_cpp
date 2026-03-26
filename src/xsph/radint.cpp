// Radial integration for multipole matrix elements.
// Converted from src/XSPH/radint.f

#include "radint.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::math {
    void csomm(const double dr[], const FeffComplex dp[], const FeffComplex dq[],
               double dpas, FeffComplex& da, int m, int np);
    double cwig3j(int j1, int j2, int j3, int m1, int m2, int m3);
}

// C++ implementation of xmult (converted from src/XSPH/xmult.f)
// Computes relativistic angular multipliers from Grant, Adv. Phys. 19, 747 (1970)
namespace {

// 6j symbol for special case j2 = j1+1 (Messiah eq. C.38, C.39)
// All input angular momenta are multiplied by 2
double sixj(int j1, int j2, int j3, int j4, int j5) {
    double aa = 0.0;
    if (j2 == j1 + 1) {
        if (j4 == j3 + 1) {
            // eq. C.38
            int g2 = j5 - 1;
            if (g2 >= std::abs(j1 - j3) && g2 <= j1 + j3) {
                aa = (1.0 + (g2 + j1 - j3) / 2.0) * (1.0 + (g2 - j1 + j3) / 2.0) /
                     static_cast<double>((j1 + 1) * (j1 + 2) * (j3 + 1) * (j3 + 2));
                aa = std::sqrt(aa);
                int sign_exp = 1 + (g2 + j1 + j3) / 2;
                if (sign_exp % 2 != 0) aa = -aa;
            }
        } else if (j3 == j4 + 1) {
            // eq. C.39
            int g2 = j5;
            if (g2 >= std::abs(j1 - j4) && g2 <= j1 + j4) {
                aa = (1.0 - (g2 - j1 - j4) / 2.0) * (2.0 + (g2 + j1 + j4) / 2.0) /
                     static_cast<double>((j1 + 1) * (j1 + 2) * (j4 + 1) * (j4 + 2));
                aa = std::sqrt(aa);
                int sign_exp = 1 + (g2 + j1 + j4) / 2;
                if (sign_exp % 2 != 0) aa = -aa;
            }
        }
    }
    return aa;
}

// 9j symbol calculation using Messiah eq. C.41
void ninej(int lam, int lamp, int ls, int j2, int jp2, int lb, double& aa) {
    if (ls > lb) {
        aa = -(ls + lb + 1) * sixj(1, 2, 2*lb, ls+lb, 2*ls) *
             sixj(2*lb, ls+lb, 2*lamp, jp2, j2) *
             sixj(ls+lb, 2*ls, 2*lam, j2, 2*lamp);
    } else if (ls < lb) {
        aa = -(ls + lb + 1) * sixj(1, 2, 2*lb, ls+lb, 2*ls) *
             sixj(ls+lb, 2*lb, jp2, 2*lamp, j2) *
             sixj(2*ls, ls+lb, j2, 2*lam, 2*lamp);
    } else {
        // ls == lb (magnetic dipole)
        aa = -(2*ls + 2) * sixj(1, 2, 2*lb, 2*lb+1, 2*lb) *
             sixj(2*lb, 2*lb+1, 2*lamp, jp2, j2) *
             sixj(2*lb, 2*lb+1, j2, 2*lam, 2*lamp);
        aa += -(2*ls) * sixj(1, 2, 2*lb, 2*lb-1, 2*lb) *
              sixj(2*lb-1, 2*lb, jp2, 2*lamp, j2) *
              sixj(2*lb-1, 2*lb, 2*lam, j2, 2*lamp);
    }
}

// xmult: compute relativistic angular multipliers (Grant eq. 6.30)
// k = ikap (final kappa), kp = kinit (initial kappa)
// ls = l (Bessel function order), lb = mult (multipole order)
// xm1, xm2 = output angular multipliers
void xmult_cpp(int k, int kp, int ls, int lb,
               feff::FeffComplex& xm1, feff::FeffComplex& xm2) {
    using feff::coni;

    // Set the factor in front of Bessel function (eq. 6.26)
    feff::FeffComplex alslb(0.0, 0.0);
    if (ls + 1 == lb) {
        // e.g. dipole and quadrupole transition
        double aa = (2*lb - 1) * (lb + 1) / 2.0;
        feff::FeffComplex coni_pow = std::pow(coni, ls);
        alslb = coni_pow * std::sqrt(aa);
    } else if (ls - 1 == lb) {
        // e.g. cross dipole-octupole
        double aa = (2*lb + 3) * lb / 2.0;
        feff::FeffComplex coni_pow = std::pow(coni, ls);
        alslb = coni_pow * std::sqrt(aa);
    } else if (ls == lb) {
        // e.g. magnetic dipole
        feff::FeffComplex coni_pow = std::pow(coni, ls);
        alslb = coni_pow * static_cast<double>(2*lb + 1) / std::sqrt(2.0);
    } else {
        alslb = 0.0;
    }

    // Set all angular momenta
    int j2 = 2 * std::abs(k) - 1;
    int a = 1;
    if (k > 0) a = -1;
    int jp2 = 2 * std::abs(kp) - 1;
    int ap = 1;
    if (kp > 0) ap = -1;

    // Calculate xm1 (beta=1 in eq. 6.30)
    int lam = (j2 - a) / 2;
    int lamp = (jp2 + ap) / 2;
    if (2*lam == j2 - a && 2*lamp == jp2 + ap) {
        double aa;
        ninej(lam, lamp, ls, j2, jp2, lb, aa);
        double sign = (lam % 2 == 0) ? 1.0 : -1.0;
        xm1 = alslb * aa * feff::math::cwig3j(lam, ls, lamp, 0, 0, 1) * sign
             * std::sqrt(6.0 * (j2+1) * (jp2+1) * (2*lb+1) * (2*lam+1) * (2*lamp+1));
        xm1 = xm1 * coni;
    } else {
        xm1 = 0.0;
    }

    // Calculate xm2 (beta=-1 in eq. 6.30)
    lam = (j2 + a) / 2;
    lamp = (jp2 - ap) / 2;
    if (2*lam == j2 + a && 2*lamp == jp2 - ap) {
        double aa;
        ninej(lam, lamp, ls, j2, jp2, lb, aa);
        double sign = (lam % 2 == 0) ? 1.0 : -1.0;
        xm2 = alslb * aa * feff::math::cwig3j(lam, ls, lamp, 0, 0, 1) * sign
             * std::sqrt(6.0 * (j2+1) * (jp2+1) * (2*lb+1) * (2*lam+1) * (2*lamp+1));
        // factor -1 due to complex conjugation of i*Q_k
        xm2 = -coni * xm2;
    } else {
        xm2 = 0.0;
    }
}

} // anonymous namespace

namespace feff::xsph {

void xrci(int mult, const FeffComplex xm[4],
          double dgc0, double dpc0, FeffComplex p, FeffComplex q,
          const double bf[3], FeffComplex& value) {
    if (mult == 0) {
        // Electric dipole with j0 and j2 contributions
        value = dgc0 * q * (xm[1] * bf[0] + xm[3] * bf[2]) +
                dpc0 * p * (xm[0] * bf[0] + xm[2] * bf[2]);
    } else {
        value = (xm[1] * dgc0 * q + xm[0] * dpc0 * p) * bf[1];
    }
}

void radint(int ifl, int mult, const double bf[][nrptx],
            int kinit, const double dgc0[], const double dpc0[],
            int ikap, FeffComplex p[], FeffComplex q[],
            FeffComplex pn[], FeffComplex qn[],
            const double ri[], double dx, int ilast, int iold,
            FeffComplex xrc[], FeffComplex xnc[],
            FeffComplex xrcold[], FeffComplex xncold[],
            FeffComplex& xirf) {

    FeffComplex temp(1.0, 0.0);
    int linit = kinit;
    if (kinit < 0) linit = -kinit - 1;
    int lfin = ikap;
    if (ikap < 0) lfin = -ikap - 1;

    // Set multipliers from Grant, Adv. Phys. 19, 747 (1970)
    FeffComplex xm[4] = {0.0, 0.0, 0.0, 0.0};
    if (ifl < 0) {
        // Nonrelativistic case
        int ji2 = 2 * std::abs(kinit) - 1;
        int jf2 = 2 * std::abs(ikap) - 1;
        if (mult == 0 || mult == 2) {
            int ll = (mult == 2) ? 2 : 1;
            int ll2 = 2 * ll;
            temp = std::sqrt(static_cast<double>((ji2 + 1) * (jf2 + 1))) *
                   feff::math::cwig3j(jf2, ll2, ji2, 1, 0, 2);
            int sign = (std::abs(ikap) % 2 == 0) ? 1 : -1;
            temp *= static_cast<double>(sign);
            int ls = ll - 1;
            FeffComplex coni_pow = std::pow(coni, ls);
            xm[0] = temp * static_cast<double>(ll2 + 1) * coni_pow * static_cast<double>(2 * ls + 1) *
                feff::math::cwig3j(ls, 1, ll, 0, 0, 1) * feff::math::cwig3j(ls, 1, ll, 0, 1, 1);
            xm[2] = 0.0;
        }
    } else if (mult == 0) {
        // Relativistic electric dipole
        xmult_cpp(ikap, kinit, 0, 1, xm[0], xm[1]);
        xmult_cpp(ikap, kinit, 2, 1, xm[2], xm[3]);
    } else {
        xm[2] = 0.0;
        xm[3] = 0.0;
        if (mult == 2) {
            xmult_cpp(ikap, kinit, 1, 2, xm[0], xm[1]);
        } else {
            xmult_cpp(ikap, kinit, 1, 1, xm[0], xm[1]);
        }
    }

    // Radial integrals
    int ia = std::abs(ifl);
    int is = ifl / ia;

    if (ia == 1) {
        // Single radial integral for rkk
        for (int i = 0; i < ilast; i++) {
            xnc[i] = 0.0;
            if (is > 0) {
                double bf_local[3] = {bf[0][i], bf[1][i], bf[2][i]};
                xrci(mult, xm, dgc0[i], dpc0[i], p[i], q[i], bf_local, xrc[i]);
            } else {
                // Nonrelativistic
                if (mult == 0) {
                    temp = xm[0] * bf[0][i] + xm[2] * bf[2][i];
                } else if (mult == 2) {
                    temp = xm[0] * bf[1][i];
                }
                temp = temp * coni;
                xrc[i] = ri[i] * (dgc0[i] * p[i] + dpc0[i] * q[i]) * temp;
            }
            if (iold == 1) xrcold[i] = xrc[i];
        }
        xirf = FeffComplex(lfin + linit + 2, 0.0);
        if (mult > 0) xirf = xirf + 1.0;
        feff::math::csomm(ri, xrc, xnc, dx, xirf, 0, ilast);
    } else {
        // Double radial integral
        if (ia == 2) {
            for (int i = 0; i < ilast; i++) {
                if (is > 0) {
                    double bf_local[3] = {bf[0][i], bf[1][i], bf[2][i]};
                    xrci(mult, xm, dgc0[i], dpc0[i], pn[i], qn[i], bf_local, xnc[i]);
                    xrci(mult, xm, dgc0[i], dpc0[i], p[i], q[i], bf_local, xrc[i]);
                } else {
                    if (mult == 0) temp = (xm[0] * bf[0][i] + xm[2] * bf[2][i]) * coni;
                    else if (mult == 2) temp = xm[0] * bf[1][i] * coni;
                    xrc[i] = ri[i] * (dgc0[i] * p[i] + dpc0[i] * q[i]) * temp;
                    xnc[i] = ri[i] * (dgc0[i] * pn[i] + dpc0[i] * qn[i]) * temp;
                }
                if (iold == 1) xncold[i] = xnc[i];
            }
        } else if (ifl == 3 && iold == 2) {
            for (int i = 0; i < ilast; i++) {
                xrc[i] = xrcold[i];
                double bf_local[3] = {bf[0][i], bf[1][i], bf[2][i]};
                xrci(mult, xm, dgc0[i], dpc0[i], pn[i], qn[i], bf_local, xnc[i]);
            }
        } else if (ifl == 4 && iold == 2) {
            for (int i = 0; i < ilast; i++) {
                double bf_local[3] = {bf[0][i], bf[1][i], bf[2][i]};
                xrci(mult, xm, dgc0[i], dpc0[i], p[i], q[i], bf_local, xrc[i]);
                xnc[i] = xncold[i];
            }
        }

        // Same for all double integrals
        if ((iold == 0 && ia == 2) || (ifl > 2 && iold == 2)) {
            int lpwr = lfin + linit + 2;
            xirf = 2.0 * xrc[0] * ri[0] / static_cast<double>(lpwr + 1);
            xnc[0] = xnc[0] * xirf;
            for (int i = 1; i < ilast; i++) {
                xirf = xirf + (xrc[i - 1] + xrc[i]) * (ri[i] - ri[i - 1]);
                xnc[i] = xnc[i] * xirf;
            }
            for (int i = 0; i < ilast; i++) {
                xrc[i] = 0.0;
            }
            xirf = FeffComplex(lpwr + 1 + linit + 1 - lfin, 0.0);
            feff::math::csomm(ri, xrc, xnc, dx, xirf, 0, ilast);
        }
    }
}

} // namespace feff::xsph
