// Sort atoms, build geometry, and create polarization tensor.
// Converted from: src/RDINP/ffsort.f
// Original coded by A.L. Ankudinov, 1998.

#include "ffsort.hpp"
#include "mkptz.hpp"
#include "../common/logging.hpp"
#include "../math/distance.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <vector>

namespace feff::rdinp {

void ffsort(int iabs, bool doptz, FeffInput& inp) {
    auto& log = feff::common::logger();

    constexpr double big = 1.0e5;

    // Local geometry arrays
    int nat = 0;
    std::array<int, nphx + 1> iatph{};
    std::vector<int> iphat(natx, 0);
    std::vector<int> index(natx, 0);
    // rat stored as [natx][3]
    std::vector<std::array<double, 3>> rat(natx, {0.0, 0.0, 0.0});

    int nabs = 1; // local copy for json output

    // The Fortran version reads from json files; in C++ the data is already
    // in inp (ratx, iphatx, natt, iphabs, rclabs, ipol, ispin, le2, etc.)

    // Find the first absorber (iphabs type) or the iabs-th atom of that type
    int iatabs = 0;
    int icount = 0;
    int ifound = 0;
    for (int iat = 0; iat < inp.natt; ++iat) {
        if (inp.iphatx[iat] == 0) inp.iphatx[iat] = inp.iphabs;
        if (inp.iphatx[iat] == inp.iphabs) icount++;
        if (ifound == 0 && icount > 0 &&
            (icount == iabs || (iabs <= 0 && icount == 1))) {
            iatabs = iat + 1; // 1-based
            ifound = 1;
        }
    }

    // Sanity checks
    if (iatabs == 0 && inp.natt > 1) {
        log.wlog(" No absorbing atom (unique pot 0 or iphabs in CFAVERAGE card) was defined.");
        throw std::runtime_error("RDINP: no absorbing atom found");
    }
    if (inp.iphabs == 0 && icount > 1) {
        log.wlog(" More than one absorbing atom (potential 0)");
        log.wlog(" Only one absorbing atom allowed");
        throw std::runtime_error("RDINP: multiple absorbers");
    }
    if ((icount > 0 && icount < nabs) || nabs <= 0) {
        nabs = icount;
        log.wlog(" Averaging over ALL atoms of iphabs type");
    }

    // Make absorbing atom first in the short list
    if (iatabs != 0) {
        rat[0] = {0.0, 0.0, 0.0};
        iphat[0] = 0;
        index[0] = iatabs;
    }

    // Build smaller list of atoms within rclabs of absorber
    nat = 1;
    for (int iat = 0; iat < inp.natt; ++iat) {
        if (iat != (iatabs - 1)) {
            // Distance between this atom and absorber
            double r0[3] = {inp.ratx[iat][0], inp.ratx[iat][1], inp.ratx[iat][2]};
            double r1[3] = {inp.ratx[iatabs - 1][0], inp.ratx[iatabs - 1][1], inp.ratx[iatabs - 1][2]};
            double tmp = feff::math::dist(r0, r1);
            if (tmp > 0.1 && tmp <= inp.rclabs) {
                nat++;
                if (nat > natx) {
                    std::ostringstream ss;
                    ss << " Number of atoms " << nat << " exceeds max allowed"
                       << " for the pathfinder = " << natx;
                    log.wlog(ss.str());
                    log.wlog(" Use or reduce rclabs in CFAVERAGE card");
                    log.wlog(" Or increase parameter natx and recompile");
                    throw std::runtime_error("RDINP: too many atoms for pathfinder");
                }
                int idx = nat - 1;
                rat[idx][0] = inp.ratx[iat][0] - inp.ratx[iatabs - 1][0];
                rat[idx][1] = inp.ratx[iat][1] - inp.ratx[iatabs - 1][1];
                rat[idx][2] = inp.ratx[iat][2] - inp.ratx[iatabs - 1][2];
                iphat[idx] = inp.iphatx[iat];
                index[idx] = iat + 1; // 1-based
            }
        }
    }

    // Sort atoms by distance from absorber (selection sort)
    for (int iat = 0; iat < nat - 1; ++iat) {
        double r2min = rat[iat][0] * rat[iat][0] + rat[iat][1] * rat[iat][1] + rat[iat][2] * rat[iat][2];
        int imin = iat;
        for (int i = iat + 1; i < nat; ++i) {
            double r2 = rat[i][0] * rat[i][0] + rat[i][1] * rat[i][1] + rat[i][2] * rat[i][2];
            if (r2 < r2min) {
                r2min = r2;
                imin = i;
            }
        }
        if (imin != iat) {
            // Swap atoms iat and imin
            std::swap(rat[iat], rat[imin]);
            std::swap(iphat[iat], iphat[imin]);
            std::swap(index[iat], index[imin]);
        }
    }

    // Rotate coordinates and make polarization tensor
    if (doptz) {
        // Convert to C-style 2D array for mkptz
        // mkptz expects double rat[][3]
        mkptz(inp.ipol, inp.elpty, inp.evec.data(), inp.xivec.data(),
              inp.ispin, inp.spvec.data(), nat,
              reinterpret_cast<double(*)[3]>(rat.data()),
              inp.angks, inp.le2,
              reinterpret_cast<FeffComplex(*)[3]>(inp.ptz.data()));
    }

    // Find model atoms for unique pots closest to absorber
    for (int iph = 1; iph <= nphx; ++iph) {
        iatph[iph] = 0;
    }
    // Absorbing atom is first in list
    iatph[0] = 1;
    inp.nph = 0;
    for (int iph = 1; iph <= nphx; ++iph) {
        double rabs = big;
        for (int iat = 1; iat < nat; ++iat) {
            if (iph == iphat[iat]) {
                double r0[3] = {rat[iat][0], rat[iat][1], rat[iat][2]};
                double r1[3] = {rat[0][0], rat[0][1], rat[0][2]};
                double tmp = feff::math::dist(r0, r1);
                if (tmp < rabs) {
                    rabs = tmp;
                    iatph[iph] = iat + 1; // 1-based
                }
            }
        }
        if (iatph[iph] > 0) inp.nph = iph;
    }

    // Check for atoms closer than 1.75 bohr (~0.93 Angstrom)
    for (int iat = 0; iat < nat; ++iat) {
        for (int jat = iat + 1; jat < nat; ++jat) {
            double r0[3] = {rat[iat][0], rat[iat][1], rat[iat][2]};
            double r1[3] = {rat[jat][0], rat[jat][1], rat[jat][2]};
            double rtmp = feff::math::dist(r0, r1);
            if (rtmp < 1.75 * bohr) {
                log.wlog(" WARNING:  TWO ATOMS VERY CLOSE TOGETHER.  CHECK INPUT.");
                int iatx = index[iat];
                int jatx = index[jat];
                {
                    std::ostringstream ss;
                    ss << " atoms " << iatx << " " << jatx;
                    log.wlog(ss.str());
                }
                {
                    std::ostringstream ss;
                    ss << std::setw(5) << iatx << std::scientific << std::setprecision(5)
                       << std::setw(13) << inp.ratx[iatx - 1][0]
                       << std::setw(13) << inp.ratx[iatx - 1][1]
                       << std::setw(13) << inp.ratx[iatx - 1][2];
                    log.wlog(ss.str());
                }
                {
                    std::ostringstream ss;
                    ss << std::setw(5) << jatx << std::scientific << std::setprecision(5)
                       << std::setw(13) << inp.ratx[jatx - 1][0]
                       << std::setw(13) << inp.ratx[jatx - 1][1]
                       << std::setw(13) << inp.ratx[jatx - 1][2];
                    log.wlog(ss.str());
                }
                log.wlog(" Run continues in case you really meant it.");
            }
        }
    }

    // Check that absorber was found
    if (iatabs <= 0) {
        log.wlog(" Absorbing atom coords not specified.");
        log.wlog(" Cannot find multiple scattering paths.");
        throw std::runtime_error("RDINP: absorbing atom not found");
    }

    // Note: The Fortran version writes geom.dat/geom.json and global.json here.
    // In the C++ version, the geometry is kept in memory. The caller can
    // write output files as needed via separate routines.
}

} // namespace feff::rdinp
