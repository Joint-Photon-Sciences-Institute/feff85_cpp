// Simplified geometric preparation without Debye-Waller factors.
// Converted from: yprep.f
// yprep is the same as xprep for negative idwopt.
// Simplifies calls in SCF and LDOS where DW factors should not enter.

#include "yprep.hpp"
#include "xprep.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

namespace feff::fms {

void yprep(int iph0, int nat, int& inclus,
           const int* iphat, float rmax, const float* rat,
           FMSData& data) {

    auto& cluster = data.cluster;
    auto& rot     = data.rotation;
    auto& lnlm    = data.lnlm;
    auto& dw      = data.dw;

    // Initialize geometrical arrays
    std::fill(cluster.xphi.begin(), cluster.xphi.end(), 0.0f);
    std::fill(cluster.xrat.begin(), cluster.xrat.end(), 0.0f);
    std::fill(cluster.iphx.begin(), cluster.iphx.end(), 0);
    inclus = 0;

    // Working copies
    std::vector<float> rat2(nat * 3);
    std::vector<int> iphat2(nat);

    // Find the central atom (ipot == iph0)
    int icen = -1;
    for (int i = 0; i < nat; ++i) {
        iphat2[i] = iphat[i];
        if (iphat[i] == iph0) {
            if (icen < 0) {
                icen = i;
            } else if (iph0 == 0) {
                throw std::runtime_error(
                    "yprep: More than one atom with ipot=0 in extended cluster");
            }
        }
    }
    if (icen < 0) {
        throw std::runtime_error("yprep: No central atom found with ipot=" +
                                 std::to_string(iph0));
    }

    // Center on central atom
    for (int i = 0; i < nat; ++i) {
        rat2[i * 3 + 0] = rat[i * 3 + 0] - rat[icen * 3 + 0];
        rat2[i * 3 + 1] = rat[i * 3 + 1] - rat[icen * 3 + 1];
        rat2[i * 3 + 2] = rat[i * 3 + 2] - rat[icen * 3 + 2];
    }

    // Sort by distance
    std::vector<double> ra(nat);
    atheap(nat, rat2.data(), iphat2.data(), ra.data());

    // Define cluster within rmax
    float rmax2 = rmax * rmax;
    inclus = 0;
    for (int i = 0; i < nat; ++i) {
        float rr = rat2[i * 3 + 0] * rat2[i * 3 + 0]
                 + rat2[i * 3 + 1] * rat2[i * 3 + 1]
                 + rat2[i * 3 + 2] * rat2[i * 3 + 2];
        if (rr > rmax2) {
            inclus = i;
            break;
        }
    }
    if (inclus == 0) inclus = nat;
    if (inclus > nclusx) inclus = nclusx;

    // Copy into cluster
    for (int iat = 0; iat < inclus; ++iat) {
        cluster.iphx[iat] = iphat2[iat];
        cluster.rat(iat, 0) = rat2[iat * 3 + 0];
        cluster.rat(iat, 1) = rat2[iat * 3 + 1];
        cluster.rat(iat, 2) = rat2[iat * 3 + 2];
    }

    // Calculate rotation matrices and phi angles
    rotint(rot);
    int lplus1 = lx + 1;
    int mplus1 = lx + 1;

    for (int i = 0; i < inclus; ++i) {
        for (int j = 0; j < inclus; ++j) {
            float xbeta, xphi_val;
            getang(cluster.xrat.data(), nclusx, i, j, xbeta, xphi_val);
            cluster.phi(i, j) = xphi_val;

            if (i != j) {
                for (int kk = 0; kk <= 1; ++kk) {
                    float beta = (kk == 1) ? -xbeta : xbeta;
                    rotxan(lplus1, mplus1, beta, i, j, kk, cluster, rot);
                }
            }
        }
    }

    // Calculate normalization factors
    xanlm(lplus1, mplus1, lnlm);

    // Zero DW factors (no thermal disorder in yprep)
    std::fill(dw.sigsqr.begin(), dw.sigsqr.end(), 0.0f);
}

} // namespace feff::fms
