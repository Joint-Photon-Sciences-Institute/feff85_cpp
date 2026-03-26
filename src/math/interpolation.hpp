#pragma once

// Polynomial interpolation and extrapolation
// Converted from src/MATH/terp.f, terpc.f, polint.f

#include <feff/types.hpp>

namespace feff::math {

// Binary search: find index of grid point immediately below x.
// xx must be monotonically increasing, size n.
// Returns 0 if x < xx[0], n-1 if x >= xx[n-1].
int locat(double x, int n, const double xx[]);

// Polynomial interpolation of order m (max 3).
// Returns interpolated y0 at x0 from arrays x[n], y[n].
void terp(const double x[], const double y[], int n, int m, double x0, double& y0);

// Complex polynomial interpolation of order m (max 3).
void terpc(const double x[], const FeffComplex y[], int n, int m,
           double x0, FeffComplex& y0);

// Polynomial interpolation through n points (Neville's algorithm).
// Returns y = P(x) and dy = error estimate.
void polint(const double xa[], const double ya[], int n, double x,
            double& y, double& dy);

// Complex version of polint.
void polinc(const double xa[], const FeffComplex ya[], int n, double x,
            FeffComplex& y, FeffComplex& dy);

// Linear interpolation for mixed single/double precision (terp1 from ff2chi).
// x and y are single precision arrays.
void terp1(const float x[], const float y[], int n, double x0, double& y0);

// Binary search for single precision array.
int locat1(double x, int n, const float xx[]);

} // namespace feff::math
