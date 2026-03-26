// Norman radius from electron density.
// Converted from src/POT/frnrm.f

#include "frnrm.hpp"
#include "../common/radial_grid.hpp"
#include "../common/logging.hpp"
#include "../math/sommerfeld.hpp"
#include <cmath>

namespace feff::pot {

void frnrm(const double* rho, int iz, double& rnrm)
{
    // rho is 1-based in Fortran; here 0-based. rho[i] => Fortran rho(i+1)
    // rr(j) is 1-based Fortran; use feff::common::rr(j) which is 1-based.
    using feff::common::rr;

    // Extended Simpson integration (see somm2) initial terms
    // Fortran indices 1..nrptx map to C++ 0..nrptx-1
    // Use r*r*r instead of std::pow(r,3) to match Fortran's integer
    // exponentiation (x**3 = x*x*x), avoiding exp(3*log(x)) rounding.
    double r1 = rr(1), r2 = rr(2), r3 = rr(3);
    double sum = (9.0 * rho[0] * r1*r1*r1 +
                  28.0 * rho[1] * r2*r2*r2 +
                  23.0 * rho[2] * r3*r3*r3) / 480.0;

    // Add initial point (r=0) correction (see somm2)
    double dpas = 0.05;
    double d1 = 3.0;
    double dd = std::exp(dpas) - 1.0;
    double db = d1 * (d1 + 1.0) * dd * std::exp((d1 - 1.0) * dpas);
    db = rr(1) / db;
    dd = rr(1) * (1.0 + 1.0 / (dd * (d1 + 1.0))) / d1;
    sum = sum + dd * rho[0] * rr(1) * rr(1) - db * rho[1] * rr(2) * rr(2);

    double r4 = rr(4), r5 = rr(5), r6 = rr(6);
    double fl = rho[3] * r4*r4*r4;
    double fr = rho[4] * r5*r5*r5;
    double frr = rho[5] * r6*r6*r6;
    sum = sum + (25.0 * fl + 12.0 * fr - frr) / 480.0;

    int inrm = 0;
    double x = 0.0;
    double sumsav = 0.0;

    // Fortran loop: do i = 7, nrptx => C++ i = 6..nrptx-1
    for (int i = 6; i < nrptx; ++i) {
        double fll = fl;
        fl = fr;
        fr = frr;
        double ri1 = rr(i + 1);
        frr = rho[i] * ri1*ri1*ri1;
        sumsav = sum;
        sum = sum + (13.0 * (fr + fl) - fll - frr) / 480.0;
        if (sum >= static_cast<double>(iz)) {
            inrm = i - 1;  // Fortran: inrm = i-2, but i is 1-based there
                            // In Fortran i starts at 7, inrm = i-2 => grid index i-2
                            // C++ i starts at 6 (=Fortran 7), inrm = (i+1)-2 = i-1
            x = (static_cast<double>(iz) - sumsav) / (sum - sumsav);
            goto found;
        }
    }

    feff::common::logger().wlog(" FRNRM Could not integrate enough charge to reach required z.");
    throw std::runtime_error("FRNRM-1");

found:
    // inrm is 0-based C++ index, but rr() takes 1-based
    // Fortran: rnrm = rr(inrm)*(1 + x*0.05)
    // inrm in Fortran is 1-based grid index; our inrm is that same value minus 1
    // Actually let's be careful: Fortran inrm = i-2 where i is 1-based (7..nrptx)
    // In C++ i is 0-based (6..nrptx-1), so Fortran_i = i+1, Fortran_inrm = i+1-2 = i-1
    // So rr(Fortran_inrm) = rr(i-1) where i-1 is still 1-based since i >= 6 => i-1 >= 5
    // Our inrm = i-1 which is 1-based Fortran index. Good.
    rnrm = rr(inrm) * (1.0 + x * 0.05);

    // Next order correction (ALA 3/97)
    double dx05 = 0.05;
    double x0 = 8.8;
    int jnrm = static_cast<int>((std::log(rnrm) + x0) / dx05) + 2;
    int i0 = jnrm + 1;

    // Build ri and xpc arrays (1-based Fortran convention)
    // We need indices 1..jnrm+2 => C++ 0..jnrm+1
    double ri[251];
    double xpc[251];
    for (int ir = 0; ir < jnrm + 2 && ir < 251; ++ir) {
        ri[ir] = rr(ir + 1);
        xpc[ir] = rho[ir] * ri[ir] * ri[ir];
    }

    double xirf = 2.0;
    feff::math::somm2(ri, xpc, dx05, xirf, rnrm, 0, i0);

    // dq is how many new electrons are within Norman sphere
    double dn1 = xirf - static_cast<double>(iz);
    // inrm here is 1-based Fortran index
    // xpc[inrm-1] = xpc(inrm) in Fortran, xpc[inrm] = xpc(inrm+1)
    double x2 = x - dn1 / ((1.0 - x) * xpc[inrm - 1] + x * xpc[inrm]);

    if (std::abs(x2 - x) > 0.0001) {
        xirf = 2.0;
        rnrm = rr(inrm) * (1.0 + x2 * 0.05);
        feff::math::somm2(ri, xpc, dx05, xirf, rnrm, 0, i0);
        double dn2 = xirf - static_cast<double>(iz);
        // Newton-Raphson method to find zeroes
        x = x2 - dn2 * (x2 - x) / (dn2 - dn1);
    }

    rnrm = rr(inrm) * (1.0 + x * 0.05);
}

} // namespace feff::pot
