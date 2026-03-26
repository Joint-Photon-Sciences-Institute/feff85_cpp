// Calculate photoelectron orbital using LDA in the Dirac equation.
// Converted from: src/FOVRG/wfirdc.f

#include "wfirdc.hpp"
#include "nucdec.hpp"
#include "potdvp.hpp"
#include "solout.hpp"
#include "solin.hpp"
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::fovrg {

void wfirdc(FeffComplex eph, int kap[], int nmax[], const FeffComplex vxc[],
            FeffComplex ps[], FeffComplex qs[], FeffComplex aps[10],
            FeffComplex aqs[10], int irr, int ic3, FeffComplex vm[],
            int jri, int& iwkb, FovrgState& state)
{
    auto& orb = state.orb;
    auto& scf = state.scf;
    auto& work = state.work;
    auto& nuclear = state.nuclear;
    auto& mesh = state.mesh;
    auto& error = state.error;
    int norb = scf.norb;
    int nz = scf.nz;
    int idim = mesh.idim;
    int ibgp = orb.ibgp;

    double cl = work.cl;  // use speed of light set by dfovrg (= alpinv)
    double dz = static_cast<double>(nz);
    work.dz = dz;

    // Make r-mesh and calculate nuclear potential
    // Match Fortran: dr1 = nz * exp(-8.8), where -8.8 is single-precision.
    // In Fortran, exp(single) returns single, integer*single = single,
    // then promoted to double on assignment to dr1.
    double dr1 = static_cast<double>(nz * std::exp(-8.8f));
    nucdec(nuclear.anoy, mesh.dr, nuclear.dvn, dz, mesh.hx,
           nuclear.nuc, idim, 10, dr1);

    // Calculate fl and fix for each orbital
    double a = (dz / cl) * (dz / cl);
    if (nuclear.nuc > 1) a = 0.0;

    for (int j = 0; j < norb; j++) {
        double b = kap[j] * kap[j] - a;
        if (j == norb - 1) b = b + (kap[j] + 1) * ic3;
        orb.fl[j] = std::sqrt(b);
        orb.fix[j] = std::pow(mesh.dr[0], orb.fl[j] - std::abs(kap[j]));
    }

    // If irregular solution, flip sign of fl and fix for last orbital
    if (irr > 0) {
        orb.fl[norb - 1] = -orb.fl[norb - 1];
        orb.fix[norb - 1] = 1.0 / orb.fix[norb - 1];
    }

    // Use LDA potential to calculate initial wavefunction
    // Fortran: do i=1,jri-1 => dv(i) = vxc(i)/cl
    // C++ 0-based: dv[0..jri-2] = vxc[0..jri-2]/cl
    for (int i = 0; i < jri - 1; i++) {
        work.dv[i] = vxc[i] / cl;
    }
    // Fortran: do i=jri,idim => dv(i) = vxc(jri+1)/cl
    // C++ 0-based: dv[jri-1..idim-1] = vxc[jri]/cl
    for (int i = jri - 1; i < idim; i++) {
        work.dv[i] = vxc[jri] / cl;
    }

    if (error.numerr != 0) return;

    // Zero exchange potentials
    for (int i = 0; i < idim; i++) {
        work.eg[i] = FeffComplex(0.0, 0.0);
        work.ep[i] = FeffComplex(0.0, 0.0);
    }
    for (int i = 0; i < ibgp; i++) {
        work.ceg[i] = FeffComplex(0.0, 0.0);
        work.cep[i] = FeffComplex(0.0, 0.0);
    }

    // Calculate potential development coefficients
    potdvp(state);
    // av(2) += (vxc(nuc) - dvn(nuc)) / cl
    // Fortran 1-based: vxc(nuc), dvn(nuc) => 0-based: vxc[nuc-1], dvn[nuc-1]
    work.av[1] = work.av[1] + (vxc[nuclear.nuc - 1] - FeffComplex(nuclear.dvn[nuclear.nuc - 1], 0.0)) / cl;

    // Set initial values for resolution of Dirac equation
    if (irr < 0) {
        if (a > 0.0) {
            aps[0] = FeffComplex(1.0, 0.0);
            if (kap[norb - 1] < 0) {
                aqs[0] = aps[0] * dz / (cl * (kap[norb - 1] - orb.fl[norb - 1]));
            } else {
                aqs[0] = aps[0] * cl * (kap[norb - 1] + orb.fl[norb - 1]) / dz;
            }
        } else {
            if (kap[norb - 1] < 0) {
                aps[0] = FeffComplex(1.0, 0.0);
                aqs[0] = FeffComplex(0.0, 0.0);
            } else {
                aps[0] = FeffComplex(0.0, 0.0);
                aqs[0] = FeffComplex(1.0, 0.0);
            }
        }
    }

    // np = 1 + int((8.8 + log(10.0))/hx)
    mesh.np = 1 + static_cast<int>((8.8 + std::log(10.0)) / mesh.hx);
    if (idim < mesh.np) mesh.np = idim;
    if (nmax[norb - 1] > mesh.np) nmax[norb - 1] = mesh.np;

    if (irr < 0) {
        solout(eph, FeffComplex(orb.fl[norb - 1], 0.0), aps[0], aqs[0],
               kap[norb - 1], jri, nmax[norb - 1], ic3, vm, iwkb,
               work, mesh);
    } else {
        solin(eph, FeffComplex(orb.fl[norb - 1], 0.0), kap[norb - 1],
              jri, nmax[norb - 1], ic3, vm, iwkb, work, mesh);
    }

    // Copy results back
    for (int i = 0; i < 10; i++) {
        aps[i] = work.ag[i];
        aqs[i] = work.ap[i];
    }
    for (int i = 0; i < idim; i++) {
        ps[i] = work.gg[i];
        qs[i] = work.gp[i];
    }
}

} // namespace feff::fovrg
