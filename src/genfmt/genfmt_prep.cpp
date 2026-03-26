// GENFMT initialization (reads phase.pad, sets up arrays).
// Converted from GENFMT/genfmt_prep.f

#include "genfmt_prep.hpp"
#include "snlm.hpp"
#include "../common/physics_utils.hpp"
#include "../common/periodic_table.hpp"
#include "../common/file_io.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <complex>
#include <string>

namespace feff::genfmt {

void genfmt_prep(const std::string& phpad, int ispin,
                 PhaseData& pd, SpinPhaseData& spd, NlmData& nlm,
                 GenfmtPrepResult& prep) {

    // ---- Read phase.pad via rdxsph (C++ equivalent) ----
    common::PhaseData xd = common::read_xsph(phpad);

    // Copy scalars into prep and pd
    prep.ne     = xd.ne;
    prep.ne1    = xd.ne1;
    prep.ne3    = xd.ne3;
    prep.npot   = xd.nph;
    prep.ihole  = xd.ihole;
    prep.rnrmav = xd.rnrmav;
    prep.xmu    = xd.xmu;
    prep.edge   = xd.edge;
    prep.ik0    = xd.ik0;
    prep.ixc    = xd.ixc;
    prep.rs     = xd.rs;
    prep.vint   = xd.vint;
    prep.lmaxp1 = xd.lmaxp1;

    pd.edge    = xd.edge;
    pd.rnrmav  = xd.rnrmav;
    pd.xmu     = xd.xmu;
    pd.ihole   = xd.ihole;

    const int ne  = xd.ne;
    const int nph = xd.nph;
    const int nsp_file = xd.nsp;  // spin channels stored in phase.pad
    const int ltot2p1 = 2 * ltot + 1;

    // Copy energy mesh
    for (int ie = 0; ie < ne; ++ie) {
        pd.em[ie] = xd.em[ie];
    }

    // Copy iz and potlbl
    for (int iph = 0; iph <= nph; ++iph) {
        pd.iz[iph] = xd.iz[iph];
        pd.potlbl[iph] = xd.potlbl[iph];
    }

    // Copy lmax(ie, iph)
    for (int iph = 0; iph <= nph; ++iph) {
        for (int ie = 0; ie < ne; ++ie) {
            pd.lmax[ie][iph] = xd.lmax[iph * ne + ie];
        }
    }

    // Copy eref2(ie, isp)
    for (int isp = 0; isp < nsp_file; ++isp) {
        for (int ie = 0; ie < ne; ++ie) {
            spd.eref2[ie][isp] = xd.eref[isp * ne + ie];
        }
    }

    // Copy ph4(ie, il, isp, iph)
    // xd.ph layout: ph[iph * nsp * ltot2p1 * ne + isp * ltot2p1 * ne + (ll+ltot) * ne + ie]
    for (int iph = 0; iph <= nph; ++iph) {
        for (int isp = 0; isp < nsp_file; ++isp) {
            for (int il = -xd.lmaxp1 + 1; il < xd.lmaxp1; ++il) {
                for (int ie = 0; ie < ne; ++ie) {
                    size_t flat = static_cast<size_t>(iph) * nsp_file * ltot2p1 * ne
                                + static_cast<size_t>(isp) * ltot2p1 * ne
                                + static_cast<size_t>(il + ltot) * ne
                                + ie;
                    ph4_at(spd, ie, il, isp, iph) = xd.ph[flat];
                }
            }
        }
    }

    // Copy rkk2(ie, kdif, isp)
    // xd.rkk layout: rkk[isp * 8 * ne + kdif * ne + ie]
    for (int isp = 0; isp < nsp_file; ++isp) {
        for (int kdif = 0; kdif < 8; ++kdif) {
            for (int ie = 0; ie < ne; ++ie) {
                spd.rkk2[ie][kdif][isp] = xd.rkk[isp * 8 * ne + kdif * ne + ie];
            }
        }
    }

    // ---- Set initial state kappa and angular momentum ----
    common::setkap(prep.ihole, prep.kinit, prep.linit);
    prep.ilinit = prep.linit + 1;

    // Number of spin channels
    prep.nsp = 1;
    if (ispin == 1) prep.nsp = nspx;

    if (prep.nsp == 1) {
        // For ispin=2 the variables are already written into is=0 positions
        int is = 0;
        for (int ie = 0; ie < prep.ne; ++ie) {
            pd.eref[ie] = spd.eref2[ie][is];
        }
        for (int iph = 0; iph <= prep.npot; ++iph) {
            for (int ie = 0; ie < prep.ne; ++ie) {
                for (int il = -pd.lmax[ie][iph]; il <= pd.lmax[ie][iph]; ++il) {
                    ph_at(pd, ie, il, iph) = ph4_at(spd, ie, il, is, iph);
                }
            }
        }
    } else {
        // Average over two spin directions
        for (int ie = 0; ie < prep.ne; ++ie) {
            pd.eref[ie] = (spd.eref2[ie][0] + spd.eref2[ie][prep.nsp - 1]) / 2.0;
        }
        for (int iph = 0; iph <= prep.npot; ++iph) {
            for (int ie = 0; ie < prep.ne; ++ie) {
                for (int il = -pd.lmax[ie][iph]; il <= pd.lmax[ie][iph]; ++il) {
                    ph_at(pd, ie, il, iph) =
                        (ph4_at(spd, ie, il, 0, iph) +
                         ph4_at(spd, ie, il, prep.nsp - 1, iph)) / 2.0;
                }
            }
        }
    }

    // Set nlm factors for use later
    snlm(ltot + 1, mtot + 1, nlm.xnlm);

    // Set potential labels from atomic symbols where blank
    // (In Fortran: atsym(iz(i)), here we skip symbol lookup
    //  since potlbl is typically already filled by rdxsph)

    // Make xk and ck arrays for later use
    for (int ie = 0; ie < prep.ne; ++ie) {
        // Real momentum (k)
        prep.xk[ie] = common::getxk(pd.em[ie].real() - prep.edge);
        // Complex momentum (p)
        prep.ck[ie] = std::sqrt(2.0 * (pd.em[ie] - pd.eref[ie]));
        prep.ckmag[ie] = std::abs(prep.ck[ie]);
        prep.xkr[ie] = prep.xk[ie];  // Real part of xk (already real)
    }

    // Central atom phase shift index
    prep.ll = prep.linit + 1;
    if (prep.kinit < 0) prep.ll = -prep.ll;

    // Initialize path counters
    prep.npath = 0;
    prep.ntotal = 0;
    prep.nused = 0;
    prep.xportx = -1.0;

}

} // namespace feff::genfmt
