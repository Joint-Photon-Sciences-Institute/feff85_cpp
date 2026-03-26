// Main pathfinder: finds multiple scattering paths using heap algorithm.
// Converted from: src/PATH/paths.f
//
// This is single precision, units are Angstroms. BE CAREFUL!

#include "paths.hpp"
#include "heap.hpp"
#include "ipack.hpp"
#include "sortix.hpp"
#include "ccrit.hpp"
#include "../common/logging.hpp"
#include "../math/distance.hpp"
#include "../par/parallel.hpp"
#include <feff/dimensions.hpp>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

namespace feff::path {

void paths(const float ckspc[], const float fbetac[], const float xlamc[],
           float pcritk, float pcrith, float critpw,
           int nncrit, float& rmax, int nlegxx, float rfms,
           int& nat, const double ratdp[][3], const int iphat[], const int ibounc[]) {

    auto& log = feff::common::logger();

    // Shared atom data -- transform to single precision, 0-based indexing
    AtomData atoms;
    atoms.resize(nat + 1);  // 0 to nat-1 after conversion (nat decremented below)

    // read_geom_json fills arrays 0-based: ratdp[0..nat-1], iphat[0..nat-1],
    // ibounc[0..nat-1].  Convert to internal atoms structure (also 0-based)
    // with nat decremented by 1 (matching the Fortran convention where the
    // original 1-based nat becomes 0-based after "nat = nat - 1").
    for (int iat = 0; iat < nat; ++iat) {
        int j = iat;
        for (int i = 0; i < 3; ++i) {
            atoms.rat[j][i] = static_cast<float>(ratdp[iat][i]);
        }
        atoms.ipot[j] = iphat[iat];
        atoms.i1b[j] = ibounc[iat];
    }
    nat = nat - 1;  // Switch to 0-based count (0..nat inclusive)
    atoms.resize(nat + 1);

    // nlegxx -> npatxx
    int npatxx = std::min(npatx, nlegxx - 1);
    // Input rmax is one-way; double it for total path length
    rmax = rmax * 2;

    int iat0 = -1;
    atoms.i1b[0] = 0;

    // Find index for the central atom (ipot==0)
    for (int iat = 0; iat <= nat; ++iat) {
        if (atoms.ipot[iat] == 0 && iat0 < 0) iat0 = iat;
    }

    // iclus: 0 inside rfms, 1 outside
    std::vector<int> iclus(nat + 1, 0);
    for (int iat = 0; iat <= nat; ++iat) {
        float rtmp = feff::math::sdist(atoms.rat[iat].data(), atoms.rat[iat0].data());
        iclus[iat] = 0;
        if (rtmp > rfms) iclus[iat] = 1;
    }

    // If central atom is not atom 0, permute
    if (iat0 != 0) {
        std::swap(atoms.rat[0], atoms.rat[iat0]);
        std::swap(atoms.ipot[0], atoms.ipot[iat0]);
        std::swap(iclus[0], iclus[iat0]);
    }

    // Count first-bounce atoms
    int n1b = 0;
    for (int i = 1; i <= nat; ++i) {
        if (atoms.i1b[i] > 0) ++n1b;
    }
    if (n1b < 1) feff::par::par_stop("At least one 1st bounce atoms required.");
    if (rmax >= big) feff::par::par_stop("Hey, get real with rmax!");

    // Make title
    char title[81];
    std::snprintf(title, sizeof(title),
        "PATH  Rmax=%6.3f,  Keep_limit=%5.2f, Heap_limit%5.2f  Pwcrit=%5.2f%%",
        rmax / 2, pcritk, pcrith, critpw);

    char slog[512];
    std::snprintf(slog, sizeof(slog), "    Rmax%8.4f  keep and heap limits%12.7f%12.7f",
        rmax / 2, pcritk, pcrith);
    log.wlog(slog);
    log.wlog("    Preparing neighbor table");

    // Build neighbor table m(-1:nat, 0:nat)
    // m[i][j] = atom index giving j-th shortest path from i back to central atom
    // i = -1 is first-bounce row, i = 0..nat are regular rows
    // Store with offset: m_idx = (i+1), so row -1 -> 0, row 0 -> 1, etc.
    int m_rows = nat + 2;  // -1 to nat
    int m_cols = nat + 1;  // 0 to nat

    // Use vector for m table
    std::vector<int> m(m_rows * m_cols, 0);
    auto m_at = [&](int i, int j) -> int& { return m[(i + 1) * m_cols + j]; };

    // Temporary arrays for sorting (1-based)
    std::vector<float> r_sort(nat + 2);    // r[1..nat+1]
    std::vector<int> mindex(nat + 2);      // mindex[1..nat+1]

    for (int i = -1; i <= nat; ++i) {
        int ir = (i == -1) ? 0 : i;
        for (int j = 0; j <= nat; ++j) {
            // r starts at element 1
            r_sort[j + 1] = feff::math::sdist(atoms.rat[ir].data(), atoms.rat[j].data());
            r_sort[j + 1] += feff::math::sdist(atoms.rat[j].data(), atoms.rat[0].data());
            // m(i,i) not needed; set to big so it sorts to end
            if (j == ir) r_sort[j + 1] = big;
            // First bounce: only allowed first bounce paths
            if (i == -1 && atoms.i1b[j] <= 0) r_sort[j + 1] = big;
        }

        // Sort row
        sortir(nat + 1, mindex.data(), r_sort.data());
        for (int j = 0; j <= nat; ++j) {
            m_at(i, j) = mindex[j + 1] - 1;
        }
    }

    // Allocate heap data structures using vectors
    // All heap arrays are 1-based (element 0 unused)
    std::vector<int> index_heap(nx + 1);
    std::vector<float> r_heap(nx + 1);
    std::vector<int> mi_heap(nx + 1);
    std::vector<int> mj_heap(nx + 1);
    std::vector<int> npat_heap(nx + 1);
    // ipat_heap[ileg][ip] -- npatx legs, nx entries
    // Flatten: ipat_heap[ip * npatx + ileg]
    std::vector<int> ipat_heap(static_cast<size_t>(npatx) * (nx + 1), 0);
    auto ipat_at = [&](int ileg, int ip) -> int& {
        return ipat_heap[ip * npatx + ileg];
    };
    std::vector<bool> keep1(nx + 1, false);

    // Initialize heap data space "next" pointers
    for (int i = 1; i < nx; ++i) {
        npat_heap[i] = i + 1;
    }
    npat_heap[nx] = -1;

    // Initial condition: first path
    int n = 1;       // Number in heap
    int nna = 0;     // Number skipped
    int nhx = n;     // Max heap size seen
    int nwrote = 0;  // Number written
    index_heap[n] = 1;
    int ip = index_heap[n];
    int next = 2;
    mi_heap[ip] = -1;
    mj_heap[ip] = 0;
    npat_heap[ip] = 1;
    ipat_at(npat_heap[ip] - 1, 1) = m_at(mi_heap[ip], mj_heap[ip]);

    // Initialize keep criterion
    float xcalcx = -1.0f;
    bool keep_flag, keep1_flag;
    float r_tmp;
    ccrit(npat_heap[ip], &ipat_at(0, ip), ckspc, fbetac, xlamc,
          rmax, pcrith, pcritk, nncrit, atoms.ipot.data(),
          r_heap[n], keep_flag, keep1_flag, xcalcx, iclus.data(), atoms);
    keep1[ip] = keep1_flag;

    // Open paths.bin for output
    std::ofstream fbin("paths.bin", std::ios::binary);
    if (!fbin) feff::par::par_stop("Cannot open paths.bin");

    // Write header
    int nhead = 1;
    fbin.write(reinterpret_cast<const char*>(&nhead), sizeof(int));
    fbin.write(title, 80);
    fbin.write(reinterpret_cast<const char*>(&nat), sizeof(int));
    for (int i = 0; i <= nat; ++i) {
        fbin.write(reinterpret_cast<const char*>(atoms.rat[i].data()), 3 * sizeof(float));
        fbin.write(reinterpret_cast<const char*>(&atoms.ipot[i]), sizeof(int));
        fbin.write(reinterpret_cast<const char*>(&atoms.i1b[i]), sizeof(int));
    }

    int np = 0;    // Paths found and saved
    int nbx = 0;   // Max npat seen
    bool ok = false;
    bool wlabel = false;

    // Main loop: "while not done"
    for (;;) {
        if (r_heap[1] > rmax || np >= npx || n <= 0) {
            if (n <= 0) ok = true;
            break;
        }

        // Save element at top of heap
        ip = index_heap[1];
        int npat0 = npat_heap[ip];
        int ipat0[8];
        for (int i = 0; i < npat0; ++i) {
            ipat0[i] = ipat_at(i, ip);
        }
        float r0 = r_heap[1];

        // Write path if last atom is not central and meets pcritk
        if (ipat0[npat0 - 1] != 0 && keep1[ip]) {
            ++np;
            int iout[3];
            ipack(iout, npat0, ipat0);
            fbin.write(reinterpret_cast<const char*>(&r0), sizeof(float));
            fbin.write(reinterpret_cast<const char*>(iout), 3 * sizeof(int));
            ++nwrote;

            if (np % 1000 == 0) {
                if (!wlabel) {
                    log.wlog("    nfound  heapsize  maxheap  maxscatt   reff");
                    wlabel = true;
                }
                std::snprintf(slog, sizeof(slog), "    %6d%9d%9d%7d%12.4f",
                    np, n, nhx, nbx, r0 / 2);
                log.wlog(slog);
            }
        }

        if (np >= npx) {
            std::snprintf(slog, sizeof(slog), "%15d paths found.  (np >= npx)", np);
            log.wlog(slog);
            break;
        }

        // Replace last atom in path from top of heap
        mj_heap[ip] = mj_heap[ip] + 1;
        if (mi_heap[ip] == -1 && atoms.i1b[m_at(mi_heap[ip], mj_heap[ip])] <= 0) {
            r_heap[1] = big;
            keep_flag = false;
        } else if (mj_heap[ip] >= nat) {
            r_heap[1] = big;
            keep_flag = false;
        } else {
            ipat_at(npat_heap[ip] - 1, ip) = m_at(mi_heap[ip], mj_heap[ip]);
            ccrit(npat_heap[ip], &ipat_at(0, ip), ckspc, fbetac, xlamc,
                  rmax, pcrith, pcritk, nncrit, atoms.ipot.data(),
                  r_heap[1], keep_flag, keep1_flag, xcalcx, iclus.data(), atoms);
            keep1[ip] = keep1_flag;
        }

        // Remove from heap if too long or fails criterion
        if (r_heap[1] > rmax || !keep_flag) {
            index_heap[1] = index_heap[n];
            r_heap[1] = r_heap[n];
            npat_heap[ip] = next;
            next = ip;
            --n;
            ++nna;
        }
        if (n > 0 && npat_heap[index_heap[1]] > nbx) {
            nbx = npat_heap[index_heap[1]];
        }

        // Bubble down
        if (n > 0) hdown(r_heap.data(), index_heap.data(), n);

        // Add an atom onto the saved path, put at end of heap
        if (npat0 + 1 <= npatxx) {
            ip = next;
            if (ip < 0) {
                // Heap full
                break;
            }
            int next0 = npat_heap[ip];
            for (int i = 0; i < npat0; ++i) {
                ipat_at(i, ip) = ipat0[i];
            }
            mi_heap[ip] = ipat0[npat0 - 1];
            mj_heap[ip] = 0;
            npat_heap[ip] = npat0 + 1;
            ipat_at(npat_heap[ip] - 1, ip) = m_at(mi_heap[ip], mj_heap[ip]);

            bool kp1tmp;
            ccrit(npat_heap[ip], &ipat_at(0, ip), ckspc, fbetac, xlamc,
                  rmax, pcrith, pcritk, nncrit, atoms.ipot.data(),
                  r_tmp, keep_flag, kp1tmp, xcalcx, iclus.data(), atoms);

            if (r_tmp > rmax || !keep_flag) {
                npat_heap[ip] = next0;
                ++nna;
            } else {
                next = next0;
                ++n;
                if (n > nhx) nhx = n;
                index_heap[n] = ip;
                r_heap[n] = r_tmp;
                keep1[ip] = kp1tmp;
                if (npat_heap[index_heap[n]] > nbx) {
                    nbx = npat_heap[index_heap[n]];
                }
                hup(r_heap.data(), index_heap.data(), n);
            }
        }
    }

    if (!ok) {
        log.wlog("   Internal path finder limit exceeded -- "
                 "path list may be incomplete.");
    }
    fbin.close();

    std::snprintf(slog, sizeof(slog), "    Paths found%9d   (maxheap, maxscatt%8d%4d)",
        np, nhx, nbx);
    log.wlog(slog);

    // Restore rmax
    rmax = rmax / 2;
}

} // namespace feff::path
