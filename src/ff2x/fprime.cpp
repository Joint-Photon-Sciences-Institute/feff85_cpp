// F-prime (anomalous scattering) calculation.
// Converted from: src/FF2X/fprime.f
// Calculates f' including solid state and lifetime effects
// using algorithm in Ankudinov, Rehr DANES paper.

#include "fprime.hpp"
#include "xscorr.hpp"  // for lorenz()

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "../math/interpolation.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>

namespace feff::ff2x {

static constexpr double eps4_local = 1.0e-4;

FeffComplex funlog(int icase, double xloss, double w, double dele) {
    FeffComplex result(0.0, 0.0);

    if (icase == 1) {
        if (std::abs(dele) >= eps4_local) {
            result = coni / 2.0 / pi *
                     (std::log(FeffComplex(-xloss, dele) / w) +
                      std::log(FeffComplex(xloss, dele) / w));
        } else {
            result = coni / pi * std::log(std::abs(xloss / w));
        }
    } else if (icase == 2) {
        if (std::abs(dele) >= eps4_local) {
            result = coni / 2.0 / pi * FeffComplex(w, xloss) * (
                std::log(FeffComplex(-xloss, dele) / w) /
                    FeffComplex(w + dele, xloss) +
                std::log(FeffComplex(xloss, dele) / w) /
                    FeffComplex(w + dele, -xloss));
        } else {
            result = coni / pi * std::log(std::abs(xloss / w)) *
                     (1.0 + coni * xloss / FeffComplex(w, -xloss));
        }
    } else if (icase == 3) {
        if (std::abs(dele) >= eps4_local) {
            result = -(w + xloss) / 2.0 / pi * (
                std::log(FeffComplex(-xloss, dele) / w) /
                    FeffComplex(dele, w + xloss) +
                std::log(FeffComplex(xloss, dele) / w) /
                    FeffComplex(dele, w - xloss));
        } else {
            result = coni / pi * std::log(std::abs(xloss / w)) *
                     (1.0 + xloss / (w - xloss));
        }
    }

    return result;
}

void fpint(const FeffComplex emxs[], const FeffComplex xmu[],
           int n1, int n2, double dele, double xloss, double eps4_val,
           double efermi, FeffComplex& value) {

    // Last interval - similar to Matsubara pole
    FeffComplex z1 = emxs[n2 - 1] - efermi;  // 0-indexed
    FeffComplex z2 = emxs[n2 - 2] - efermi;
    value = -coni / pi * (z1 - dele) /
            (xloss * xloss + std::pow(z1 - dele, 2)) *
            xmu[n2 - 1] * (2.0 * (z1 - z2));

    // All other intervals (Fortran n1 to n2-2, 1-indexed -> C++ n1-1 to n2-3)
    for (int i = n1 - 1; i <= n2 - 3; ++i) {
        z1 = emxs[i] - efermi;
        z2 = emxs[i + 1] - efermi;
        FeffComplex bb = (xmu[i + 1] * (z2 - dele) - xmu[i] * (z1 - dele)) /
                         xloss / (z2 - z1);
        FeffComplex aa = xmu[i] * (z1 - dele) / xloss - bb * z1;
        FeffComplex c1 = (aa + bb * (dele + coni * xloss)) / 2.0 / coni;

        if (std::abs(dele - z1.real()) < eps4_val &&
            std::abs(dele - z2.real()) < eps4_val) {
            value += -coni / pi * c1 *
                     std::log(std::abs((z2 - dele - coni * xloss) /
                                       (z1 - dele - coni * xloss)));
        } else {
            value += -coni / pi * c1 *
                     std::log((z2 - dele - coni * xloss) /
                              (z1 - dele - coni * xloss));
        }

        c1 = -(aa + bb * (dele - coni * xloss)) / 2.0 / coni;
        value += -coni / pi * c1 *
                 std::log((z2 - dele + coni * xloss) /
                          (z1 - dele + coni * xloss));
    }
}

void fpintp(const double em[], const FeffComplex xmu[],
            int n2, double dele, double xloss, double efermi,
            FeffComplex& value) {

    value = FeffComplex(0.0, 0.0);

    // All intervals (1 to n2-1 in Fortran, 0 to n2-2 in C++)
    for (int i = 0; i < n2 - 1; ++i) {
        double x1 = em[i] - efermi;
        double x2 = em[i + 1] - efermi;
        double de = (x2 - x1) / 2.0;
        double x0 = (em[i] + em[i + 1]) / 2.0;

        FeffComplex aa;
        math::terpc(em, xmu, n2, 3, x0, aa);

        FeffComplex bb = (xmu[i + 1] - xmu[i]) / (x2 - x1);
        FeffComplex cc = (xmu[i + 1] - aa - bb * de) / (de * de);

        FeffComplex z1_val = dele - x0 + efermi - coni * xloss;
        FeffComplex z2_val = dele - x0 + efermi + coni * xloss;

        value += 2.0 * de * bb + 2.0 * z1_val * de * cc +
                 std::log((de - z1_val) / (-de - z1_val)) *
                 (aa + bb * z1_val + cc * z1_val * z1_val);
        value += 2.0 * de * bb + 2.0 * z2_val * de * cc +
                 std::log((de - z2_val) / (-de - z2_val)) *
                 (aa + bb * z2_val + cc * z2_val * z2_val);
    }

    // Tail of xmu to infinity approximated by aa/(w-bb)^2
    double x1 = em[n2 - 2];
    double x2 = em[n2 - 1];
    double a = std::sqrt(std::abs(xmu[n2 - 2].real() / xmu[n2 - 1].real()));
    double b = (a * x1 - x2) / (a - 1.0);
    if (b > x1) b = 0.0;
    FeffComplex aa_tail = xmu[n2 - 1] * std::pow(x2 - b, 2);
    FeffComplex z1_val = dele - coni * xloss - b;
    FeffComplex z2_val = dele + coni * xloss - b;
    double x0 = x2 - b;

    value += std::log(x0 / (x0 - z1_val)) * aa_tail / (z1_val * z1_val) -
             aa_tail / z1_val / x0;
    value += std::log(x0 / (x0 - z2_val)) * aa_tail / (z2_val * z2_val) -
             aa_tail / z2_val / x0;

    // Multiply by constant factor
    value = -coni / 2.0 / pi * value;
}

void fprime(double ei, FeffComplex emxs[], int ne1, int ne3, int ne, int ik0,
            FeffComplex xsec[], double xsnorm[], FeffComplex chia[],
            double vrcorr, double vicorr, FeffComplex cchi[]) {

    double efermi = emxs[ne1].real();  // ne1+1 in Fortran, ne1 in 0-indexed
    double xloss = emxs[0].imag();
    int ne2 = ne - ne1 - ne3;

    // Build xmu in complex energy plane
    FeffComplex xmu[nex];
    if (ne2 > 0) {
        // DANES
        for (int ie = 0; ie < ne1; ++ie)
            xmu[ie] = coni * static_cast<double>(xsnorm[ie]) + xsnorm[ie] * chia[ie];
        for (int ie = ne1; ie < ne1 + ne2; ++ie)
            xmu[ie] = xsnorm[ie] * chia[ie];
        for (int ie = ne - ne3; ie < ne; ++ie)
            xmu[ie] = coni * static_cast<double>(xsnorm[ie]);
    } else {
        // FPRIME
        for (int ie = 0; ie < ne; ++ie)
            xmu[ie] = xsec[ie] + xsnorm[ie] * chia[ie];
    }

    // Handle vrcorr shift
    if (std::abs(vrcorr) > eps4_local) {
        FeffComplex bb_val = xmu[ik0];
        efermi -= vrcorr;
        double omega_tmp[nex];
        for (int ie = 0; ie < ne1; ++ie)
            omega_tmp[ie] = emxs[ie].real();
        math::terpc(omega_tmp, xmu, ne1, 1, efermi, bb_val);
        for (int ie = 0; ie < ne2; ++ie)
            emxs[ne1 + ie] -= vrcorr;
        if (std::abs(xmu[ik0]) > eps4_local)
            bb_val = bb_val / xmu[ik0];
        for (int ie = ne1; ie < ne - ne3; ++ie)
            xmu[ie] = xmu[ie] * bb_val;
    }

    // Handle vicorr broadening
    if (vicorr > eps4_local) {
        xloss += vicorr;
        double omega_tmp[nex];
        for (int ie = 0; ie < ne2; ++ie)
            omega_tmp[ie] = emxs[ne1 + ie].imag();
        FeffComplex aa;
        math::terpc(omega_tmp, &xmu[ne1], ne2, 1, xloss, aa);
        for (int ie = 0; ie < ne1; ++ie) {
            double xx = vicorr * vicorr /
                        (vicorr * vicorr + std::pow(emxs[ie].real() - efermi, 2));
            xmu[ie] = xmu[ie] * (1.0 - xx) + aa * xx;
            emxs[ie] += coni * vicorr;
        }
    }

    // Output diagnostic arrays
    double dout[7][nex];

    // Cycle over energy points on horizontal grid
    for (int ie = 0; ie < ne1; ++ie) {
        dout[0][ie] = emxs[ie].real() * hart;
        double dele = emxs[ie].real() - efermi;
        double delp = -dele - 2.0 * ei;

        cchi[ie] = FeffComplex(0.0, 0.0);

        if (ne2 > 0) {
            if (std::abs(dele) < eps4_local) dele = 0.0;
            double w1 = emxs[ne1].imag();
            double w2 = emxs[ne1 + 1].imag();
            double w3 = emxs[ne1 + 2].imag();

            // Matsubara pole
            FeffComplex temp = lorenz(xloss, w1, dele) * xmu[ne1] * 2.0 * coni * w1;
            temp += lorenz(xloss, w1, delp) * xmu[ne1] * 2.0 * coni * w1;
            dout[1][ie] = temp.real();

            // Sommerfeld correction
            temp = coni * w1 * w1 / 6.0 *
                   (lorenz(xloss, w3, dele) * xmu[ne1 + 2] -
                    lorenz(xloss, w2, dele) * xmu[ne1 + 1]) / (w3 - w2);
            dout[2][ie] = temp.real();

            cchi[ie] = lorenz(xloss, w1, dele) * xmu[ne1] * 2.0 * coni * w1 +
                       coni * w1 * w1 / 6.0 *
                       (lorenz(xloss, w3, dele) * xmu[ne1 + 2] -
                        lorenz(xloss, w2, dele) * xmu[ne1 + 1]) / (w3 - w2);

            // Negative pole contribution
            cchi[ie] += lorenz(xloss, w1, delp) * xmu[ne1] * 2.0 * coni * w1 +
                        coni * w1 * w1 / 6.0 *
                        (lorenz(xloss, w3, delp) * xmu[ne1 + 2] -
                         lorenz(xloss, w2, delp) * xmu[ne1 + 1]) / (w3 - w2);

            // Theta function contribution
            if (dele < eps4_local) cchi[ie] -= xmu[ie];
            if (std::abs(dele) < eps4_local) cchi[ie] += xmu[ie] / 2.0;

            // Anomalous contribution
            temp = FeffComplex(0.0, 0.0);
            double wp_fp = 2.0 * ei;
            if (dele >= eps4_local) temp = xmu[ie];
            if (std::abs(dele) < eps4_local) temp = xmu[ie] / 2.0;
            temp = temp + xmu[ik0] * funlog(1, xloss, wp_fp, dele);
            dout[3][ie] = temp.real();

            // Integration over vertical axis
            int n1 = ne1 + 1;   // Fortran ne1+2, 0-indexed = ne1+1
            int n2 = ne - ne3;
            FeffComplex val;
            fpint(emxs, xmu, n1 + 1, n2, dele, xloss, eps4_local, efermi, val);
            cchi[ie] += val;
            fpint(emxs, xmu, n1 + 1, n2, delp, xloss, eps4_local, efermi, val);
            cchi[ie] += val;
        }

        // Integration over horizontal axis
        FeffComplex temp_h(0.0, 0.0);
        double emp[nex];
        FeffComplex xmup[nex];
        int n2_h;

        if (ne2 > 0) {
            // DANES
            int n1_d = ne1 - ik0;  // ne1 - ik0 + 1 in 1-indexed, but array start differs
            for (int i = ik0; i < ne1; ++i) {
                emp[i - ik0] = emxs[i].real();
                xmup[i - ik0] = coni * static_cast<double>(xsnorm[i]);
            }
            for (int i = 0; i < ne3; ++i) {
                emp[n1_d + i] = emxs[ne - ne3 + i].real();
                xmup[n1_d + i] = xmu[ne - ne3 + i];
            }
            n2_h = n1_d + ne3;
        } else {
            // FPRIME
            int n1_f = 0;
            for (int i = 0; i < ne1; ++i) {
                if (n1_f == 0 && emxs[i].real() > emxs[ne1].real())
                    n1_f = i;
            }
            for (int i = 0; i < ne3; ++i) {
                emp[i] = emxs[ne1 + i].real();
                xmup[i] = xmu[ne1 + i];
            }
            n2_h = ne3;
        }

        FeffComplex val_h;
        fpintp(emp, xmup, n2_h, dele, xloss, efermi, val_h);
        temp_h += val_h;
        fpintp(emp, xmup, n2_h, delp, xloss, efermi, val_h);
        temp_h += val_h;

        dout[4][ie] = temp_h.real();
        cchi[ie] += temp_h;

        // Total contribution
        FeffComplex total = xmu[ie] + cchi[ie];
        dout[5][ie] = total.real();
        dout[6][ie] = dout[5][ie] - dout[3][ie];
    }

    // Restore the input energy mesh
    if (vicorr > eps4_local) {
        for (int ie = 0; ie < ne1; ++ie)
            emxs[ie] -= coni * vicorr;
    }
    if (std::abs(vrcorr) > eps4_local) {
        for (int ie = 0; ie < ne2; ++ie)
            emxs[ne1 + ie] += vrcorr;
    }

    // Write danes.dat diagnostic file
    {
        std::ofstream fout("danes.dat");
        if (fout.is_open()) {
            fout << "# E  matsub. sommerf. anomal. tale, total, differ.\n";
            for (int ie = 0; ie < ne1; ++ie) {
                for (int j = 0; j < 7; ++j) {
                    fout << " " << std::scientific << std::setprecision(4)
                         << std::setw(12) << dout[j][ie];
                }
                fout << "\n";
            }
        }
    }
}

} // namespace feff::ff2x
