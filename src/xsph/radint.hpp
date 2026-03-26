#pragma once
// Radial integration for multipole matrix elements.
// Converted from src/XSPH/radint.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Perform radial integration for multipole matrix element or cross-section.
///
/// @param ifl    Flag: 1=rkk matrix element, -1=nonrelativistic rkk,
///               2=xsec (regular+irregular), 3=cross term irregular, 4=cross term regular
/// @param mult   Multipole: 0=E1, 1=M1, 2=E2
/// @param bf     Bessel functions for x-ray k-vector, bf[3][nrptx]
/// @param kinit  Initial state kappa
/// @param dgc0   Large Dirac component of initial orbital [nrptx]
/// @param dpc0   Small Dirac component of initial orbital [nrptx]
/// @param ikap   Final state kappa
/// @param p      Regular solution large component [nrptx]
/// @param q      Regular solution small component [nrptx]
/// @param pn     Irregular solution large component [nrptx]
/// @param qn     Irregular solution small component [nrptx]
/// @param ri     Radial grid [nrptx]
/// @param dx     Grid step
/// @param ilast  Last integration point
/// @param iold   Cross-term storage flag: 0=none, 1=store, 2=use
/// @param xrc    Work array for regular coupling [nrptx]
/// @param xnc    Work array for irregular coupling [nrptx]
/// @param xrcold Stored regular coupling [nrptx] (in/out)
/// @param xncold Stored irregular coupling [nrptx] (in/out)
/// @param xirf   Output: value of radial integral
void radint(int ifl, int mult, const double bf[][nrptx],
            int kinit, const double dgc0[], const double dpc0[],
            int ikap, FeffComplex p[], FeffComplex q[],
            FeffComplex pn[], FeffComplex qn[],
            const double ri[], double dx, int ilast, int iold,
            FeffComplex xrc[], FeffComplex xnc[],
            FeffComplex xrcold[], FeffComplex xncold[],
            FeffComplex& xirf);

/// R-dependent multipole matrix element (before r-integration).
void xrci(int mult, const FeffComplex xm[4],
          double dgc0, double dpc0, FeffComplex p, FeffComplex q,
          const double bf[3], FeffComplex& value);

} // namespace feff::xsph
