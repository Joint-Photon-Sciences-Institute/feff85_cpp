// Module 4 entry point: pathfinder.
// Converted from: src/PATH/ffmod4.f

#include "ffmod4.hpp"
#include "repath.hpp"
#include "prcrit.hpp"
#include "paths.hpp"
#include "pathsd.hpp"
#include "path_data.hpp"
#include "../common/logging.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace feff::path {

void ffmod4() {
    feff::par::par_begin();
    if (feff::par::state().worker) {
        feff::par::par_barrier();
        feff::par::par_end();
        return;
    }

    // Open log file
    auto& log = feff::common::logger();
    log.open("log4.dat");

    // Read input
    int ms, mpath, ipr4, nncrit, nlegxx, ipol, ispin;
    float pcritk, pcrith, rmax, rfms2, critpw;
    int nat;
    double rat[natx + 1][3];  // 1-based
    int iphat[natx + 1];
    int ibounc[natx + 1];
    double evec[3], xivec[3];

    repath(ms, mpath, ipr4, pcritk, pcrith, nncrit, rmax,
           nlegxx, rfms2, critpw,
           nat, rat, iphat, ibounc,
           ipol, ispin, evec, xivec);

    int eels = 0;
    if (nspx > 1) ispin = std::abs(ispin);

    if (ms == 1 && mpath == 1) {
        log.wlog(" Preparing plane wave scattering amplitudes...");

        // Allocate criterion arrays
        int fbetac_size = fbetac_dim1 * fbetac_dim2 * necrit;
        int fbeta_size  = fbetac_dim1 * fbetac_dim2 * nex;
        std::vector<float> fbetac(fbetac_size, 0.0f);
        std::vector<float> fbeta(fbeta_size, 0.0f);
        std::vector<float> ckspc(necrit, 0.0f);
        std::vector<float> xlamc(necrit, 0.0f);
        std::vector<float> cksp(nex, 0.0f);
        std::vector<float> xlam(nex, 0.0f);
        std::string potlbl[nphx + 1];

        int ne, ik0;
        prcrit(ne, nncrit, ik0,
               cksp.data(), fbeta.data(),
               ckspc.data(), fbetac.data(),
               potlbl, xlam.data(), xlamc.data());

        // Debug output (ipr4 >= 3)
        if (ipr4 >= 3 && ipr4 != 5) {
            for (int iph = 0; iph <= 1; ++iph) {
                for (int ie = 0; ie < nncrit; ++ie) {
                    char fname[32];
                    std::snprintf(fname, sizeof(fname), "fbeta%dp%d.dat", ie + 1, iph);
                    std::ofstream fdbg(fname);
                    if (fdbg) {
                        char buf[256];
                        std::snprintf(buf, sizeof(buf),
                            "# iph, ie, ckspc(ie) %5d%5d%20.6e\n"
                            "#  angle(degrees), fbeta/|p|,  fbeta",
                            iph, ie + 1, ckspc[ie]);
                        fdbg << buf << "\n";
                        for (int ibeta = -nbeta; ibeta <= nbeta; ++ibeta) {
                            float cosb = 0.025f * ibeta;
                            if (cosb > 1.0f)  cosb = 1.0f;
                            if (cosb < -1.0f) cosb = -1.0f;
                            float angle = std::acos(cosb);
                            float fb = fbetac[fbetac_idx(ibeta, iph, ie)];
                            std::snprintf(buf, sizeof(buf), "%10.4f%15.6e%15.6e",
                                angle * feff::raddeg, fb / ckspc[ie], fb);
                            fdbg << buf << "\n";
                        }
                    }
                }
            }
        }

        log.wlog(" Searching for paths...");
        paths(ckspc.data(), fbetac.data(), xlamc.data(),
              pcritk, pcrith, critpw, nncrit, rmax, nlegxx, rfms2,
              nat, rat, iphat, ibounc);

        log.wlog(" Eliminating path degeneracies...");
        pathsd(ckspc.data(), fbetac.data(), xlamc.data(),
               ne, ik0, cksp.data(), fbeta.data(), xlam.data(),
               critpw, ipr4, nncrit, potlbl,
               ipol, ispin, evec, xivec, eels);

        log.wlog(" Done with module 4: pathfinder.");
    }

    log.close();

    feff::par::par_barrier();
    feff::par::par_end();
}

} // namespace feff::path
