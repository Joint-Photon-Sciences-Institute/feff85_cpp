#pragma once

// Dimension parameters from FEFF
// Converted from src/HEADERS/dim.h

namespace feff {

// Max number of atoms for FMS
inline constexpr int nclusx = 100;

// Max number of atoms for tdlda module
inline constexpr int nclxtd = 100;

// Max number of spins: 1 for spin average; 2 for spin-dep
inline constexpr int nspx = 1;

// Max number of atoms in problem for pathfinder
inline constexpr int natx = 1000;

// Max number of atoms in problem for rdinp and ffsort
inline constexpr int nattx = 1000;

// Max orbital momentum for FMS module
inline constexpr int lx = 4;

// Max number of unique potentials (nphx must be ODD)
inline constexpr int nphx = 11;

// Max number of angular momentum (arrays 1:ltot+1)
inline constexpr int ltot = 24;

// Loucks r grid used through overlap and in phase work arrays
inline constexpr int nrptx = 1251;

// Number of energy points genfmt, etc.
inline constexpr int nex = 150;

// Max number of distinct lambda values for genfmt
inline constexpr int lamtot = 15;

// Vary mmax and nmax independently
inline constexpr int mtot = 4;
inline constexpr int ntot = 2;

// Max number of path atoms, used in path finder
inline constexpr int npatx = 8;

// Matches path finder, used in GENFMT
inline constexpr int legtot = npatx + 1;

// Max number of overlap shells (OVERLAP card)
inline constexpr int novrx = 8;

// Max number of header lines
inline constexpr int nheadx = 30;

// Max number of poles for HL multipole self energy
inline constexpr int MxPole = 1000;

} // namespace feff
