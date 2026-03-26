#pragma once
// Parallel SCF iteration stub.
// Converted from: src/POT/scmtmp.f (MPI version)
//
// Currently just delegates to the serial scmt() since we are
// running in sequential mode. Full MPI parallelization can be
// added later.

#include <feff/dimensions.hpp>
#include <complex>

namespace feff::pot {

/// Parallel SCF iteration (stub -- delegates to serial scmt).
///
/// The additional parameters (gtr, xrhole, xrhoce, yrhole, yrhoce)
/// are work arrays that would be distributed across MPI ranks in
/// the parallel version.
void scmtmp(bool verbose, int npr, int iscmt, double ecv, int nph, int nat,
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
            std::complex<float>* gtr,
            std::complex<double>* xrhole,
            std::complex<double>* xrhoce,
            std::complex<double>* yrhole,
            std::complex<double>* yrhoce);

} // namespace feff::pot
