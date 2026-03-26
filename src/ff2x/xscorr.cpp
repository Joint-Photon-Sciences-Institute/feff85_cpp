// Cross-section correction via contour integration in the complex energy plane.
// Converted from: src/FF2X/xscorr.f

#include "xscorr.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "../math/interpolation.hpp"

#include <cmath>
#include <complex>

namespace feff::ff2x {

static constexpr double eps4_local = 1.0e-4;

FeffComplex lorenz(double xloss, double w, double dele) {
    return xloss / pi / (xloss * xloss + std::pow(coni * w - dele, 2));
}

double astep(double xloss, double dele) {
    double val = 0.5 + std::atan(dele / xloss) / pi;
    if (val < 0.0) val = 0.0;
    if (val > 1.0) val = 1.0;
    return val;
}

void xscorr(int ispec, FeffComplex emxs[], int ne1, int ne, int ik0,
            FeffComplex xsec[], const double xsnorm[], FeffComplex chia[],
            double vrcorr, double vicorr, FeffComplex cchi[]) {

    int ne2 = ne - ne1;
    double efermi = emxs[ne - 1].real();
    double xloss = emxs[0].imag();
    vicorr = 0;

    // xmu = analytic function in complex energy plane
    FeffComplex xmu[nex];
    for (int ie = 0; ie < ne; ++ie) {
        xmu[ie] = xsec[ie] + xsnorm[ie] * chia[ie];
    }

    // Real frequencies
    double omega[nex];
    for (int ie = 0; ie < ne1; ++ie) {
        omega[ie] = emxs[ie].real();
    }

    FeffComplex bb(1.0, 0.0);
    if (std::abs(vrcorr) > eps4_local) {
        // Account for the Fermi level shift
        bb = xmu[ik0];
        efermi = efermi - vrcorr;
        math::terpc(omega, xmu, ne1, 1, efermi, bb);

        // Shift the vertical axis
        for (int ie = 0; ie < ne2; ++ie) {
            emxs[ne1 + ie] -= vrcorr;
        }

        // Rescale values on vertical axis
        bb = bb / xmu[ik0];
        for (int ie = ne1; ie < ne; ++ie) {
            xmu[ie] = xmu[ie] * bb;
        }
    }

    // Construct the integration contour C
    FeffComplex ec[nex], fc[nex];
    int nc = 0;

    // Points on vertical axis below xloss
    for (int ie = 0; ie < ne2; ++ie) {
        if (emxs[ne1 + ie].imag() < xloss) {
            ec[nc] = emxs[ne1 + ie];
            fc[nc] = xmu[ne1 + ie];
            nc++;
        }
    }

    // Corner at efermi + i*xloss
    int ic0 = nc;
    if (std::abs(vrcorr) > eps4_local) {
        ec[nc] = FeffComplex(efermi, xloss);
        fc[nc] = bb * xmu[ik0];
    } else {
        ec[nc] = emxs[ik0];
        fc[nc] = xmu[ik0];
    }
    nc++;

    // Points on horizontal axis above efermi (or below for emission)
    if (ispec != 2) {
        for (int ie = 0; ie < ne1; ++ie) {
            if (emxs[ie].real() - efermi > eps4_local) {
                ec[nc] = emxs[ie];
                fc[nc] = xmu[ie];
                nc++;
            }
        }
    } else {
        // Emission: points below E_fermi
        for (int ie = ne1 - 1; ie >= 0; --ie) {
            if (efermi - emxs[ie].real() > eps4_local) {
                ec[nc] = emxs[ie];
                fc[nc] = xmu[ie];
                nc++;
            }
        }
    }

    // Cycle over frequency points
    FeffComplex ff[nex];
    for (int ie = 0; ie < ne1; ++ie) {
        FeffComplex xmu0;
        if (omega[ie] >= efermi) {
            xmu0 = xmu[ie];
            if (ispec == 2) xmu0 = xmu[ik0] * bb;
        } else {
            xmu0 = xmu[ik0] * bb;
            if (ispec == 2) xmu0 = xmu[ie];
        }

        FeffComplex e1 = FeffComplex(omega[ie], xloss);
        FeffComplex e2 = FeffComplex(omega[ie], -xloss);

        for (int ic = 0; ic < nc; ++ic) {
            ff[ic] = fc[ic] - xmu0;
        }

        double dele = omega[ie] - efermi;
        cchi[ie] = xmu0 * astep(xloss, dele);
        if (ispec == 2) cchi[ie] = xmu0 - cchi[ie];

        FeffComplex corr(0.0, 0.0);

        if (std::abs(dele) < eps4_local) dele = 0.0;
        double w1 = ec[0].imag();

        // Half Matsubara pole contribution
        corr += lorenz(xloss, w1, dele) * ff[0] * coni * w1;

        // Cycle over contour points
        for (int ic = 0; ic < nc - 1; ++ic) {
            FeffComplex z1 = ec[ic];
            FeffComplex z2 = ec[ic + 1];
            FeffComplex f1 = ff[ic];
            FeffComplex f2 = ff[ic + 1];

            // Correction from pole above real axis
            FeffComplex aa(0.0, 0.0);
            if (std::abs(z1 - e1) > eps4_local && std::abs(z2 - e1) > eps4_local) {
                aa = std::log((z2 - e1) / (z1 - e1)) *
                     (f1 * (z2 - e1) + f2 * (e1 - z1));
            }
            // Second pole
            aa -= std::log((z2 - e2) / (z1 - e2)) *
                  (f1 * (z2 - e2) + f2 * (e2 - z1));

            corr += aa / (z2 - z1) / 2.0 / pi / coni;
        }

        if (ispec == 2) corr = -corr;

        cchi[ie] += corr;
        // Return result of convolution minus bare value
        cchi[ie] -= xmu[ie];
    }

    // Restore the input energy mesh
    if (std::abs(vrcorr) > eps4_local) {
        for (int ie = ne1; ie < ne; ++ie) {
            emxs[ie] += vrcorr;
        }
    }
}

} // namespace feff::ff2x
