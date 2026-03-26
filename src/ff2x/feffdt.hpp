#pragma once
// Write feffNNNN.dat and files.dat output files.
// Converted from: src/FF2X/feffdt.f, fdthea.f, fdtarr.f, fdtxdi.f

#include "ff2x_types.hpp"

#include <complex>
#include <string>
#include <vector>

namespace feff::ff2x {

/// Write feffNNNN.dat files and files.dat for all paths in the list.
/// Replaces Fortran subroutine feffdt.
void feffdt(int ntotal, const int iplst[], int nptot,
            int ntext, const std::string text[],
            int ne, int iorder, int ilinit, float rnrmav, float edge,
            const std::string potlbl[], const int iz[],
            const std::complex<float> phc[], const std::complex<float> ck[],
            const float xk[], const int index[], const int nleg[],
            const float deg[], const float reff[], const float crit[],
            const std::vector<std::vector<int>>& ipot,
            const std::vector<std::vector<std::array<float, 3>>>& rat,
            const std::vector<std::vector<float>>& achi,
            const std::vector<std::vector<float>>& phchi);

/// Build header lines for feffNNNN.dat output.
/// Replaces Fortran subroutine fdthea.
void fdthea(int ntext, const std::string text[], int ip, int iorder,
            int nleg_val, float deg, float reff, float rnrmav, float edge,
            const std::array<float, 3> rat_leg[], const int ipot_leg[],
            const int iz[], const std::string potlbl[],
            int& nlines, std::vector<std::string>& lines);

/// Compute the 7 columns of feffNNNN.dat data.
/// Replaces Fortran subroutine fdtarr.
void fdtarr(int ne, float reff, int lzero,
            const float achi[], const float phchi[],
            const std::complex<float> caps[], const float xk[],
            const std::complex<float> ck[],
            double col1[], double col2[], double col3[], double col4[],
            double col5[], double col6[], double col7[]);

} // namespace feff::ff2x
