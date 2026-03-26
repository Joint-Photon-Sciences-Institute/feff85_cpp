#pragma once
// Main SCF driver for single-configuration Dirac-Fock atom.
// Converted from: src/ATOM/scfdat.f
//
// This is the primary entry point for atomic structure calculations,
// called by the POT module to compute free-atom potentials and densities.

#include "atom_types.hpp"
#include <feff/dimensions.hpp>

namespace feff::atom {

/// Main SCF driver for Dirac-Fock atom calculation.
/// Performs iterative self-consistent field calculation to convergence.
///
/// Key outputs:
///   dgc/dpc: Dirac spinor components for all orbitals and potentials
///   adgc/adpc: development coefficients
///   srho/dmag: electron density and spin density
///   vcoul: Coulomb potential
///   eatom: total atomic energy
///   s02: S0^2 shakeup amplitude
///   eorb/kappa: orbital energies and quantum numbers
///
/// Replaces Fortran scfdat().
void scfdat(int ipr1, int iph, int nph, int iz, int ihole, double xion,
            int iunf, double vcoul[251], double srho[251], double dmag[251],
            double srhovl[251], int ispinr,
            double dgc0[251], double dpc0[251],
            double* dgc, double* dpc,     // dgc(251,30,0:nphx), flat
            double* adgc, double* adpc,   // adgc(10,30,0:nphx), flat
            double& s02, double& efrozn, double& eatom,
            double xntot[], double xnval[30],
            int indorb[], int& norbp, double eorb[30], int kappa_out[30]);

} // namespace feff::atom
