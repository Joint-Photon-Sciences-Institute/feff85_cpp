// Parallel SCF iteration stub.
// Converted from: src/POT/scmtmp.f (MPI version)
//
// Currently delegates to the serial scmt() since we are running
// in sequential mode. The extra work arrays (gtr, xrhole, etc.)
// are unused in this stub.

#include "scmtmp.hpp"
#include "scmt.hpp"

namespace feff::pot {

void scmtmp(bool verbose, int /*npr*/, int iscmt, double ecv, int nph, int nat,
            double* vclap, double* edens, double* edenvl,
            double* vtot, double* vvalgs,
            double* rmt, double* rnrm, double* qnrm,
            int ixc, double rhoint, double vint, double& xmu, int jumprm,
            double xnferm, double* xnvmu, double* xnval,
            double x0, double* ri, double dx,
            double* xnatph, double* xion, int iunf, int* iz,
            double* adgc, double* adpc, double* dgc, double* dpc,
            int ihole,
            double* rat, int* iatph, int* iphat,
            int* lmaxsc, double* rhoval, double* xnmues,
            bool& ok,
            double rgrd, int nohole, int nscmt, int icoul,
            double ca1, float rfms1, int lfms1,
            std::complex<float>* /*gtr*/,
            std::complex<double>* /*xrhole*/,
            std::complex<double>* /*xrhoce*/,
            std::complex<double>* /*yrhole*/,
            std::complex<double>* /*yrhoce*/)
{
    // Sequential stub: just call the serial version
    scmt(verbose, iscmt, ecv, nph, nat, vclap, edens,
         edenvl, vtot, vvalgs, rmt, rnrm, qnrm,
         ixc, rhoint, vint, xmu, jumprm,
         xnferm, xnvmu, xnval,
         x0, ri, dx, xnatph, xion, iunf, iz,
         adgc, adpc, dgc, dpc, ihole,
         rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok,
         rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1);
}

} // namespace feff::pot
