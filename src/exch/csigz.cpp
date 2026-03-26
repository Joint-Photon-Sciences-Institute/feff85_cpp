// Non-self-consistent self-energy with renormalization factor Z
// Converted from src/EXCH/csigz.f

#include "csigz.hpp"
#include "csigma.hpp"   // sigma1, hfexc

#include <cmath>

#include <feff/types.hpp>
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

namespace feff::exch {

void csigz(FeffComplex energy, double mu, double rs,
           double& resig, double& imsig, FeffComplex& ztot,
           const double wpscl[], const double ampfac[])
{
    constexpr int mx_iter = 1;

    double kfermi = fa / rs;
    double efermi = kfermi * kfermi / 2.0;

    FeffComplex sig_tot(0.0, 0.0);
    FeffComplex dsgde(0.0, 0.0);
    FeffComplex sigma_f(0.0, 0.0);
    double gam = 0.0;
    ztot = FeffComplex(0.0, 0.0);

    FeffComplex ckf(0.0, 0.0);
    FeffComplex ck0(0.0, 0.0);

    // Self-consistency loop (MxIter = 1)
    for (int i2 = 0; i2 < mx_iter; ++i2) {

        // Loop over poles to get SigmaF
        for (int i1 = 0; i1 < MxPole; ++i1) {
            if (wpscl[i1] < 0.0) break;

            double wp = std::sqrt(3.0 / (rs * rs * rs)) * wpscl[i1];

            ckf = FeffComplex(kfermi * 1.00001, 0.0);
            FeffComplex rel_en(efermi, 0.0);
            sigma_f += sigma1(ckf, rel_en, wp, gam, ampfac[i1], kfermi, efermi);
        }

        dsgde = FeffComplex(0.0, 0.0);

        // Loop over poles
        for (int i1 = 0; i1 < MxPole; ++i1) {
            if (wpscl[i1] < 0.0) break;

            double wp = std::sqrt(3.0 / (rs * rs * rs)) * wpscl[i1];

            // Start with ck0 = sqrt(Re(Energy) - Mu + EFermi)
            FeffComplex rel_en(energy.real() - mu + efermi, 0.0);
            ck0 = std::sqrt(2.0 * rel_en.real());

            // Find Sigma0
            FeffComplex sigma0 = sigma1(ck0, rel_en, wp, gam, ampfac[i1],
                                        kfermi, efermi);

            // Compute numerical derivative: (Sigma(E*0.001) - Sigma(E)) / (E*0.001 - E)
            FeffComplex rel_en_p = rel_en * 0.001;
            FeffComplex sigma_p = sigma1(ck0, rel_en_p, wp, gam, ampfac[i1],
                                         kfermi, efermi);

            dsgde += (sigma_p - sigma0) / (rel_en_p - rel_en);

            sig_tot += sigma0;
        }
    }

    // Add Hartree-Fock part of delta Sigma
    FeffComplex del_hf = hfexc(ck0, efermi, kfermi) - hfexc(ckf, efermi, kfermi);
    sig_tot = sig_tot - sigma_f + del_hf;

    // Form ZTot and return Re and Im parts
    ztot = 1.0 / (1.0 - dsgde);
    sig_tot = ztot * sig_tot;

    resig = sig_tot.real();
    imsig = sig_tot.imag();
}

} // namespace feff::exch
