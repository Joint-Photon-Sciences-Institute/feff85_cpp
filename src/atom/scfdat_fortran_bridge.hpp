#pragma once
// Bridge to call the original Fortran scfdat subroutine from C++ code.
// Used to test whether the C++ ATOM module reproduces the Fortran SCF results.
//
// When USE_FORTRAN_ATOM is defined, pot.cpp calls scfdat_fortran() instead
// of the C++ feff::atom::scfdat(). This links the entire Fortran ATOM module
// (scfdat.o + libfeffatom.a dependencies) to isolate whether the SCF
// convergence difference is in the ATOM module.

#include <feff/dimensions.hpp>

namespace feff::atom {

/// Call the Fortran scfdat_ subroutine via extern "C".
/// Arguments match the C++ scfdat() signature exactly.
/// Internally initializes Fortran COMMON /parallel/ before the first call.
void scfdat_fortran(
    int ipr1, int iph, int nph, int iz, int ihole, double xion,
    int iunf, double vcoul[251], double srho[251], double dmag[251],
    double srhovl[251], int ispinr,
    double dgc0[251], double dpc0[251],
    double* dgc, double* dpc,
    double* adgc, double* adpc,
    double& s02, double& efrozn, double& eatom,
    double xntot[], double xnval[30],
    int indorb[], int& norbp, double eorb[30], int kappa_out[30]);

} // namespace feff::atom
