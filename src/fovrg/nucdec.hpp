#pragma once
// Nuclear potential construction for FOVRG module.
// Converted from: src/FOVRG/nucdec.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::fovrg {

/// Construct the nuclear potential and its development coefficients at origin.
///
/// @param av    Development coefficients at origin (size 10)
/// @param dr    Radial mesh tabulation points (size nrptx)
/// @param dv    Nuclear potential output (size nrptx)
/// @param dz    Nuclear charge
/// @param hx    Exponential mesh step
/// @param nuc   Index of nuclear radius (1 = point charge) [in/out]
/// @param np    Number of tabulation points
/// @param ndor  Number of development coefficients
/// @param dr1   First tabulation point * nz [in/out]
void nucdec(double av[10], double dr[], double dv[], double dz,
            double hx, int& nuc, int np, int ndor, double& dr1);

} // namespace feff::fovrg
