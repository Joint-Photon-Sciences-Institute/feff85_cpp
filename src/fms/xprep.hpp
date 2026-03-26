#pragma once

// Geometric preparation for FMS calculation.
// Converted from: xprep.f and helper routines in xstaff.f
//   atheap  — heap-sort atoms by distance from central atom
//   getang  — polar angles between atom pairs
//   rotxan  — rotation matrix elements for atom pairs
//   rotint  — initialize rotation matrix cache
//   xanlm   — Legendre polynomial normalization factors
//   xfctst  — factorial table for normalization

#include "fms_types.hpp"

namespace feff::fms {

/// Full geometric preparation including Debye-Waller factors.
/// Replaces Fortran subroutine xprep.
///
/// @param iph0     Potential index of the central atom (0 for absorber)
/// @param idwopt   Debye-Waller option (0=CD, 1=EM, 2=RM, negative=none)
/// @param nat      Number of atoms in the extended cluster
/// @param[out] inclus  Number of atoms in the FMS cluster
/// @param npot     Number of unique potentials
/// @param iphat    Potential index for each atom in extended cluster [nat]
/// @param rmax     FMS cluster radius (Bohr)
/// @param rat      Coordinates (Bohr) [nat][3]
/// @param izx      Atomic numbers per potential [nphx+1]
/// @param rnrmav   Average Norman radius
/// @param temper   Temperature (K)
/// @param thetad   Debye temperature (K)
/// @param sig2     Global sigma^2 addend
/// @param minv     Matrix inversion method
/// @param rdirec   Direct interaction cutoff distance (Bohr)
/// @param data     FMS shared state (cluster, rotation, lnlm, dw, cg updated)
void xprep(int iph0, int idwopt, int nat, int& inclus, int npot,
           const int* iphat, float rmax, const float* rat,
           const int* izx, float rnrmav, float temper, float thetad, float sig2,
           int minv, float rdirec,
           FMSData& data);

// --- Internal helper routines (exposed for testing) ---

/// Heap-sort atoms by distance from origin.
/// Sorts rat (column-major [3][nat]) and iphat in place.
void atheap(int nat, float* rat, int* iphat, double* ra);

/// Compute polar angles (theta, phi) of vector R_i - R_j.
void getang(const float* xrat, int nclusx_dim, int i, int j,
            float& theta, float& phi);

/// Compute rotation matrix elements for atom pair (i,j).
/// k=0 forward, k=1 backward rotation.
void rotxan(int lxp1, int mxp1, float betax, int i, int j, int k,
            ClusterData& cluster, RotationData& rot);

/// Initialize rotation matrix cache.
void rotint(RotationData& rot);

/// Compute Legendre polynomial normalization factors.
void xanlm(int lmaxp1, int mmaxp1, LegendreNorm& lnlm);

/// Initialize factorial table (used by xanlm).
void xfctst(float& afac, float& flzero, float* flg);

} // namespace feff::fms
