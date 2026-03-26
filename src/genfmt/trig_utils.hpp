#pragma once

// Trigonometric utility functions for path geometry.
// Converted from GENFMT/trig.f

#include <feff/types.hpp>
#include <complex>

namespace feff::genfmt {

/// Compute cos(theta), sin(theta), cos(phi), sin(phi) for vector (x,y,z).
/// Convention: if x=y=0 and z>0, phi=0; if x=y=z=0, theta=0.
void trig(double x, double y, double z,
          double& ct, double& st, double& cp, double& sp);

/// Compute the argument (angle) of a complex number c.
/// If c is near zero, returns fi as the default angle.
void arg(FeffComplex c, double fi, double& th);

} // namespace feff::genfmt
