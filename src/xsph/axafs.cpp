// Extract AXAFS from atomic cross-section.
// Converted from src/XSPH/axafs.f

#include "axafs.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace feff::math { double determ(double* array, int nord, int nrows); }

namespace feff::xsph {

void axafs(const FeffComplex em[], double emu, const FeffComplex xsec[],
           int ne1, int ik0) {

    double xn = 0.0;
    int mm = 1;
    int np = ne1 - ik0;
    double ef = emu;

    double ee[nex], xmu_arr[nex], wt[nex];

    for (int ie = 0; ie < np; ie++) {
        ee[ie] = std::real(em[ik0 + ie] - em[ik0]) + emu;
        xmu_arr[ie] = std::imag(xsec[ik0 + ie]) * std::pow(ee[ie], xn);
    }
    for (int ie = 0; ie < np; ie++) {
        if (ie == 0) {
            wt[ie] = (ee[ie + 1] - ef) * std::pow(std::abs(ee[ie] - ef), mm);
        } else if (ie == np - 1) {
            wt[ie] = (ee[ie] - ee[ie - 1]) * std::pow(ee[ie] - ef, mm);
        } else {
            wt[ie] = (ee[ie + 1] - ee[ie - 1]) * std::pow(ee[ie] - ef, mm);
        }
    }

    double xx[5] = {0, 0, 0, 0, 0};
    double yy[3] = {0, 0, 0};
    for (int ie = 0; ie < np; ie++) {
        for (int i = 0; i < 5; i++) xx[i] += wt[ie] * std::pow(ee[ie], i);
        for (int i = 0; i < 3; i++) yy[i] += wt[ie] * xmu_arr[ie] * std::pow(ee[ie], i);
    }

    // Compute coefficients via Cramer's rule
    double xm[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            xm[i + 3 * j] = xx[i + j];
    double denom = feff::math::determ(xm, 3, 3);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            xm[i + 3 * j] = xx[i + j];
    for (int i = 0; i < 3; i++) xm[i] = yy[i];
    double aa = feff::math::determ(xm, 3, 3) / denom;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            xm[i + 3 * j] = xx[i + j];
    for (int i = 0; i < 3; i++) xm[i + 3] = yy[i];
    double bb = feff::math::determ(xm, 3, 3) / denom;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            xm[i + 3 * j] = xx[i + j];
    for (int i = 0; i < 3; i++) xm[i + 6] = yy[i];
    double cc = feff::math::determ(xm, 3, 3) / denom;

    // Normalization at edge+100 eV
    double eee = ee[0] + 100.0 / hart;
    double xnorm = (aa + bb * eee + cc * eee * eee) / std::pow(eee, xn);

    std::ofstream fout("axafs.dat");
    fout << "# File contains AXAFS. See manual for details.\n";
    fout << "#--------------------------------------------------------------\n";
    fout << "#  e, e(wrt edge), k, mu_at=(1+chi_at)*mu0_at, mu0_at, chi_at @#\n";
    fout << std::scientific << std::setprecision(5);
    for (int ie = 0; ie < np; ie++) {
        double xmu_val = std::imag(xsec[ie + ik0]);
        double xmu0 = (aa + bb * ee[ie] + cc * ee[ie] * ee[ie]) / std::pow(ee[ie], xn);
        double chiat = (xmu_val - xmu0) / xmu0;
        double eee_val = ee[ie] - ef;
        double xk;
        if (eee_val >= 0.0) {
            xk = std::sqrt(2.0 * eee_val) / bohr;
        } else {
            xk = -std::sqrt(-2.0 * eee_val) / bohr;
        }
        fout << std::fixed << std::setprecision(3) << std::setw(11) << ee[ie] * hart
             << std::setw(11) << (ee[ie] - emu) * hart
             << std::setw(8) << xk
             << std::scientific << std::setprecision(5)
             << std::setw(13) << xmu_val / xnorm
             << std::setw(13) << xmu0 / xnorm
             << std::setw(13) << chiat << "\n";
    }
}

} // namespace feff::xsph
