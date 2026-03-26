#pragma once
// Energy mesh generation for phase shift calculations.
// Converted from src/XSPH/phmesh.f
//
// Creates complex energy grid with horizontal and vertical components
// for XANES/EXAFS/DANES/FPRIME calculations.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Generate energy mesh for phase shift calculation.
///
/// @param iprint  Print level (>=3 writes emesh.dat)
/// @param ispec   Calculation type: 0=EXAFS, 1-3=XANES/XES/DANES, 4=FPRIME
/// @param edge    Chemical potential: xmu - vr0 (Hartrees)
/// @param emu     Edge energy (Hartrees)
/// @param vi0     Constant imaginary potential (Hartrees)
/// @param gamach  Core-hole broadening (Hartrees)
/// @param xkmax   Maximum k value (inverse Bohr) or emin for FPRIME
/// @param xkstep  k-step for XANES (inverse Bohr) or emax for FPRIME
/// @param vixan   Energy step for FMS calculations (Hartrees)
/// @param ne      Output: total number of energy points
/// @param ne1     Output: number of horizontal grid points
/// @param em      Output: complex energy grid [nex]
/// @param ik0     Output: grid index where k=0
/// @param ne3     Output: number of auxiliary horizontal points (DANES/FPRIME)
void phmesh(int iprint, int ispec, double edge, double emu,
            double vi0, double gamach,
            double xkmax, double xkstep, double vixan,
            int& ne, int& ne1, FeffComplex em[], int& ik0, int& ne3);

} // namespace feff::xsph
