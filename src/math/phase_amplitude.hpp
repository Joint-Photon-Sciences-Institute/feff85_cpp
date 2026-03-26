#pragma once

// Phase shift and amplitude calculations
// Converted from src/MATH/phamp.f

#include <feff/types.hpp>

namespace feff::math {

// Calculate phase shift at muffin-tin radius.
// Given wavefunction values pu,qu at rmt, momentum ck, Bessel functions
// jl,nl,jlp,nlp, and kappa sign ikap, returns phase shift ph and amplitude amp.
void phamp(double rmt, FeffComplex pu, FeffComplex qu, FeffComplex ck,
           FeffComplex jl, FeffComplex nl, FeffComplex jlp, FeffComplex nlp,
           int ikap, FeffComplex& ph, FeffComplex& amp);

// Complex arctangent: phx = atan(temp) for complex numbers
void atancc(FeffComplex temp, FeffComplex& phx);

// Complex atan2: find ampl, phx such that a = ampl*cos(phx), b = ampl*sin(phx)
void atan2c(FeffComplex a, FeffComplex b, FeffComplex& ampl, FeffComplex& phx);

} // namespace feff::math
