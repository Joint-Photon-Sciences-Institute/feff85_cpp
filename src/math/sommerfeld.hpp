#pragma once

// Simpson and extended Simpson integration on exponential grids
// Converted from src/MATH/somm.f, somm2.f, csomm.f, csomm2.f

#include <feff/types.hpp>

namespace feff::math {

// Simpson integration of (dp+dq)*dr^m from 0 to r=dr[np-1]
// dpas = exponential step; near r=0, (dp+dq) ~ cte*r^da
// Result returned in da (overwritten).
void somm(const double dr[], const double dp[], const double dq[],
          double dpas, double& da, int m, int np);

// Extended Simpson integration of dp*dr^m from 0 to r=rnrm
// with proper end corrections. dpas = exponential step.
// Near r=0, dp ~ cte*r^da. Result returned in da.
void somm2(const double dr[], const double dp[],
           double dpas, double& da, double rnrm, int m, int np);

// Complex version of somm: Simpson integration of (dp+dq)*dr^m
// dp, dq are complex arrays. da is complex (overwritten with result).
void csomm(const double dr[], const FeffComplex dp[], const FeffComplex dq[],
           double dpas, FeffComplex& da, int m, int np);

// Complex version of somm2: extended Simpson integration of dp*dr
// from 0 to r=rnrm. da is complex (overwritten with result).
void csomm2(const double dr[], const FeffComplex dp[],
            double dpas, FeffComplex& da, double rnrm, int np);

} // namespace feff::math
