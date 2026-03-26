// Bridge to call the Fortran-compiled scfdat subroutine.
//
// This file:
//   1. Declares the Fortran scfdat_ symbol as extern "C"
//   2. Declares the Fortran COMMON /parallel/ block and initializes it
//   3. Provides scfdat_fortran() which converts C++ call convention
//      (scalars by value/reference) to Fortran convention (all by pointer)
//
// The Fortran scfdat_ uses COMMON blocks internally but all I/O goes
// through its argument list, so no COMMON block marshalling is needed
// for the caller (unlike the soldir/intdir wrappers).

#include "scfdat_fortran_bridge.hpp"
#include <cstdio>

// -----------------------------------------------------------------------
// Fortran COMMON /parallel/ — must be initialized before scfdat_ is called.
// Layout from HEADERS/parallel.h:
//   common /parallel/ numprocs, my_rank, this_process,
//                     master, worker, parallel_run, par_type
//
// Fortran LOGICAL is 4 bytes on gfortran (default -fdefault-integer-4).
// The layout is:
//   numprocs   : integer (4 bytes)
//   my_rank    : integer (4 bytes)
//   this_process : integer (4 bytes)
//   master     : logical (4 bytes)
//   worker     : logical (4 bytes)
//   parallel_run : logical (4 bytes)
//   par_type   : integer (4 bytes)
// -----------------------------------------------------------------------
extern "C" {

struct FortranParallel {
    int numprocs;
    int my_rank;
    int this_process;
    int master;        // Fortran LOGICAL: 1 = .TRUE., 0 = .FALSE.
    int worker;
    int parallel_run;
    int par_type;
};

// COMMON /timing/ — also referenced by scfdat.o
struct FortranTiming {
    double wall_comm;
    double time_comm;
};

extern FortranParallel parallel_;
extern FortranTiming timing_;

// The Fortran scfdat subroutine — all arguments passed by reference
void scfdat_(
    int* ipr1, int* iph, int* nph, int* iz, int* ihole, double* xion,
    int* iunf, double* vcoul, double* srho, double* dmag, double* srhovl,
    int* ispinr, double* dgc0, double* dpc0, double* dgc, double* dpc,
    double* adgc, double* adpc,
    double* s02, double* efrozn, double* eatom,
    double* xntot, double* xnval, int* indorb, int* norbp,
    double* eorb, int* kappa);

} // extern "C"


namespace feff::atom {

// One-time initialization flag
static bool fortran_parallel_initialized = false;

static void init_fortran_parallel() {
    if (fortran_parallel_initialized) return;

    parallel_.numprocs    = 1;
    parallel_.my_rank     = 0;
    parallel_.this_process = 0;
    parallel_.master      = 1;   // .TRUE.
    parallel_.worker      = 0;   // .FALSE.
    parallel_.parallel_run = 0;  // .FALSE.
    parallel_.par_type    = 0;

    timing_.wall_comm = 0.0;
    timing_.time_comm = 0.0;

    fortran_parallel_initialized = true;
    std::fprintf(stderr, "[BRIDGE] Fortran COMMON /parallel/ initialized (master=TRUE)\n");
}

void scfdat_fortran(
    int ipr1, int iph, int nph, int iz, int ihole, double xion,
    int iunf, double vcoul[251], double srho[251], double dmag[251],
    double srhovl[251], int ispinr,
    double dgc0[251], double dpc0[251],
    double* dgc, double* dpc,
    double* adgc, double* adpc,
    double& s02, double& efrozn, double& eatom,
    double xntot[], double xnval[30],
    int indorb[], int& norbp, double eorb[30], int kappa_out[30])
{
    // Ensure Fortran runtime state is initialized
    init_fortran_parallel();

    std::fprintf(stderr, "[BRIDGE] Calling Fortran scfdat_: iph=%d iz=%d ihole=%d xion=%.4f\n",
                 iph, iz, ihole, xion);

    // Fortran convention: all arguments passed by reference (pointer).
    // Create local copies of scalar inputs that need to be passed by pointer.
    int f_ipr1 = ipr1;
    int f_iph  = iph;
    int f_nph  = nph;
    int f_iz   = iz;
    int f_ihole = ihole;
    double f_xion = xion;
    int f_iunf = iunf;
    int f_ispinr = ispinr;

    // Call Fortran scfdat_
    // Note: Fortran uses 1-based iph, and so does the C++ caller (pot.cpp
    // already passes 1-based iph values for Fortran compatibility).
    scfdat_(
        &f_ipr1, &f_iph, &f_nph, &f_iz, &f_ihole, &f_xion,
        &f_iunf, vcoul, srho, dmag, srhovl,
        &f_ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,
        &s02, &efrozn, &eatom,
        xntot, xnval, indorb, &norbp, eorb, kappa_out);

    std::fprintf(stderr, "[BRIDGE] Fortran scfdat_ returned: eatom=%.6f norbp=%d\n",
                 eatom, norbp);
}

} // namespace feff::atom
