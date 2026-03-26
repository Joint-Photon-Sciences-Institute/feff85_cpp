// Path degeneracy checker.
// Converted from: src/PATH/pathsd.f

#include "pathsd.hpp"
#include "ipack.hpp"
#include "sortix.hpp"
#include "timrep.hpp"
#include "mpprmd.hpp"
#include "outcrt.hpp"
#include "../common/logging.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

namespace feff::path {

void pathsd(const float ckspc[], const float fbetac[], const float xlamc[],
            int ne, int ik0, const float cksp[],
            const float fbeta[], const float xlam[],
            float critpw, int ipr2, int nncrit,
            const std::string potlbl[],
            int ipol, int ispin, const double evec[3], const double xivec[3],
            int eels) {

    auto& log = feff::common::logger();

    char slog[512];
    std::snprintf(slog, sizeof(slog), "    Plane wave chi amplitude filter%7.2f%%", critpw);
    log.wlog(slog);

    // Read atoms from paths.bin
    std::ifstream fbin("paths.bin", std::ios::binary);
    if (!fbin) feff::par::par_stop("Cannot open paths.bin for reading");

    int nhead;
    fbin.read(reinterpret_cast<char*>(&nhead), sizeof(int));
    std::vector<std::string> head(nhead);
    for (int i = 0; i < nhead; ++i) {
        char buf[80];
        fbin.read(buf, 80);
        head[i] = std::string(buf, 80);
    }

    int nat;
    fbin.read(reinterpret_cast<char*>(&nat), sizeof(int));

    AtomData atoms;
    atoms.resize(nat + 1);
    for (int i = 0; i <= nat; ++i) {
        fbin.read(reinterpret_cast<char*>(atoms.rat[i].data()), 3 * sizeof(float));
        fbin.read(reinterpret_cast<char*>(&atoms.ipot[i]), sizeof(int));
        fbin.read(reinterpret_cast<char*>(&atoms.i1b[i]), sizeof(int));
    }

    // Initialize counters
    int nptot = 0;
    int nuptot = 0;
    int ngs = 0;
    float xportx = eps5;
    int ndegx = -1;
    float xcalcx = -1.0f;

    // Open paths.dat output
    std::ofstream fpdat;
    if (ipr2 != 5) {
        fpdat.open("paths.dat");
        if (!fpdat) feff::par::par_stop("Cannot open paths.dat");
        for (int i = 0; i < nhead; ++i) {
            fpdat << head[i] << "\n";
        }
        for (int i = 0; i < 71; ++i) fpdat << '-';
        fpdat << "\n";
    }

    // Open crit.dat
    std::ofstream fcrit;
    if (ipr2 >= 1) {
        fcrit.open("crit.dat");
        if (!fcrit) feff::par::par_stop("Cannot open crit.dat");
        for (int i = 0; i < nhead; ++i) {
            fcrit << head[i] << "\n";
        }
        std::snprintf(slog, sizeof(slog), " Plane wave chi amplitude filter%7.2f%%", critpw);
        fcrit << slog << "\n";
        for (int i = 0; i < 71; ++i) fcrit << '-';
        fcrit << "\n";
        fcrit << " ipath nleg ndeg     r       pwcrit    "
                 "xkeep   accuracy   xheap    accuracy\n";
    }

    // Allocate path storage for one distance range
    std::vector<std::array<int, 3>> iout_arr(np1x);
    std::vector<int> index_sort(np1x + 1);
    std::vector<double> dhash(np1x + 1);  // 1-based for sortid

    // Read first path
    float r0;
    int iout0[3];
    fbin.read(reinterpret_cast<char*>(&r0), sizeof(float));
    fbin.read(reinterpret_cast<char*>(iout0), 3 * sizeof(int));
    if (fbin.eof()) goto label_999;

    // Process each total path length range
    {
    bool last = false;

    while (true) {
        ++ngs;
        float rcurr = r0;
        int np = 1;
        for (int i = 0; i < 3; ++i) iout_arr[0][i] = iout0[i];

        // Read paths with same total length
        while (true) {
            fbin.read(reinterpret_cast<char*>(&r0), sizeof(float));
            fbin.read(reinterpret_cast<char*>(iout0), 3 * sizeof(int));
            if (fbin.eof()) {
                last = true;
                break;
            }
            if (std::abs(r0 - rcurr) < eps3) {
                ++np;
                if (np > np1x) feff::par::par_stop("np > np1x");
                for (int i = 0; i < 3; ++i) iout_arr[np - 1][i] = iout0[i];
            } else {
                break;  // r0 is for next set
            }
        }

        // Hash each path
        for (int ip_i = 0; ip_i < np; ++ip_i) {
            int npat_tmp = npatx;
            int ipat_tmp[npatx];
            upack(iout_arr[ip_i].data(), npat_tmp, ipat_tmp);

            float rx[npatx], ry[npatx], rz[npatx];
            // dhash is 1-based for sortid
            timrep(npat_tmp, ipat_tmp, rx, ry, rz, dhash[ip_i + 1],
                   ipol, ispin, evec, xivec, eels, atoms);
        }

        // Sort by hash
        sortid(np, index_sort.data(), dhash.data());

        // Find ranges with same hash key
        int i0 = 1;  // 1-based index into sorted array
        while (i0 <= np) {
            int i1 = np;
            double dcurr = dhash[index_sort[i0]];
            for (int ip_i = i0 + 1; ip_i <= np; ++ip_i) {
                if (dhash[index_sort[ip_i]] != dcurr) {
                    i1 = ip_i - 1;
                    break;
                }
            }

            // Sum degeneracy and check hash collisions
            int npat0 = npatx;
            int ipat0[npatx + 1];
            float rx0[npatx], ry0[npatx], rz0[npatx];
            upack(iout_arr[index_sort[i0] - 1].data(), npat0, ipat0);
            double ddum;
            timrep(npat0, ipat0, rx0, ry0, rz0, ddum,
                   ipol, ispin, evec, xivec, eels, atoms);

            int ndeg = 0;
            for (int ii = i0; ii <= i1; ++ii) {
                int npat_tmp = npatx;
                int ipat_tmp[npatx + 1];
                float rx_tmp[npatx], ry_tmp[npatx], rz_tmp[npatx];
                upack(iout_arr[index_sort[ii] - 1].data(), npat_tmp, ipat_tmp);

                int ndpath = atoms.i1b[ipat_tmp[0]];
                timrep(npat_tmp, ipat_tmp, rx_tmp, ry_tmp, rz_tmp, ddum,
                       ipol, ispin, evec, xivec, eels, atoms);
                ndeg += ndpath;

                // Check for hash collisions
                bool ldiff = false;
                if (npat_tmp != npat0) {
                    ldiff = true;
                } else {
                    for (int iat = 0; iat < npat_tmp; ++iat) {
                        if (atoms.ipot[ipat_tmp[iat]] != atoms.ipot[ipat0[iat]]) {
                            ldiff = true;
                            break;
                        }
                    }
                    if (!ldiff) {
                        for (int ileg = 0; ileg < npat_tmp; ++ileg) {
                            if (std::abs(rx_tmp[ileg] - rx0[ileg]) > eps3 ||
                                std::abs(ry_tmp[ileg] - ry0[ileg]) > eps3 ||
                                std::abs(rz_tmp[ileg] - rz0[ileg]) > eps3) {
                                ldiff = true;
                                break;
                            }
                        }
                    }
                }
                if (ldiff) {
                    log.wlog(" WARNING!!  Two non-degenerate paths,"
                             " hashed to the same hash key!!");
                    feff::par::par_stop("hash error");
                }
            }

            // Compute path importance factors
            float xport, xheap_val, xheapr, xkeep;
            outcrt(npat0, ipat0, ckspc, nncrit, fbetac, xlamc,
                   ne, ik0, cksp, fbeta, xlam,
                   atoms.ipot.data(),
                   xport, xheap_val, xheapr, xkeep, xcalcx, atoms);

            if (xportx * ndegx <= 0) {
                xportx = xport;
                ndegx = ndeg;
            }
            float frac = 100.0f * ndeg * xport / (ndegx * xportx);

            // Write output if path is important enough
            if (frac >= critpw) {
                nptot += ndeg;
                ++nuptot;

                // Compute double-precision path parameters for output
                double rid[npatx + 1], betad[npatx + 1], etad[npatx + 1];
                mpprmd(npat0, ipat0, rid, betad, etad, atoms);

                // Write paths.dat
                if (ipr2 != 5) {
                    std::snprintf(slog, sizeof(slog),
                        " %4d%5d%8.3f  index, nleg, degeneracy, r=%8.4f",
                        nuptot, npat0 + 1, static_cast<float>(ndeg), rcurr / 2);
                    fpdat << slog << "\n";
                    fpdat << "      x           y           z     ipot  "
                             "label      rleg      beta        eta\n";
                    for (int i = 0; i < npat0; ++i) {
                        int iat = ipat0[i];
                        std::snprintf(slog, sizeof(slog),
                            "%12.6f%12.6f%12.6f%4d '%-6s' %10.4f%10.4f%10.4f",
                            atoms.rat[iat][0], atoms.rat[iat][1], atoms.rat[iat][2],
                            atoms.ipot[iat], potlbl[atoms.ipot[iat]].c_str(),
                            rid[i], betad[i] * feff::raddeg, etad[i] * feff::raddeg);
                        fpdat << slog << "\n";
                    }
                    // Central atom line
                    std::snprintf(slog, sizeof(slog),
                        "%12.6f%12.6f%12.6f%4d '%-6s' %10.4f%10.4f%10.4f",
                        atoms.rat[0][0], atoms.rat[0][1], atoms.rat[0][2],
                        atoms.ipot[0], potlbl[atoms.ipot[0]].c_str(),
                        rid[npat0], betad[npat0] * feff::raddeg,
                        etad[npat0] * feff::raddeg);
                    fpdat << slog << "\n";
                }

                // Write crit.dat
                if (ipr2 >= 1) {
                    float cmpk = xkeep * ndeg / ndegx;
                    cmpk = 100.0f - 100.0f * (std::abs(frac - cmpk) / frac);
                    float cmph;
                    if (xheap_val < 0) {
                        cmph = 100.0f;
                    } else {
                        cmph = xheap_val * ndeg / ndegx;
                        cmph = 100.0f - 100.0f * (std::abs(frac - cmph) / frac);
                    }
                    std::snprintf(slog, sizeof(slog),
                        "%6d%4d%6d%10.4f%10.4f%10.4f%8.2f%10.4f%14.3e",
                        nuptot, npat0 + 1, ndeg, rcurr / 2, frac,
                        xkeep, cmpk, xheap_val, cmph);
                    fcrit << slog << "\n";
                }
            }

            i0 = i1 + 1;
        }

        if (last) break;
    }
    }

label_999:
    if (ipr2 != 5 && fpdat.is_open()) fpdat.close();
    // Delete paths.bin
    fbin.close();
    std::remove("paths.bin");
    if (fcrit.is_open()) fcrit.close();

    std::snprintf(slog, sizeof(slog), "    Unique paths%7d,  total paths%8d",
        nuptot, nptot);
    log.wlog(slog);

    if (nuptot > 1200) {
        log.wlog(" You have found more than 1200 paths.  Genfmt");
        log.wlog(" could require a lot of time and more than 6 meg of");
        log.wlog(" storage.  Suggest a larger critpw to reduce number");
        log.wlog(" of paths.  To continue this calculation, restart");
        log.wlog(" with current paths.dat and module genfmt (3rd module");
        log.wlog(" on CONTROL card).");
        feff::par::par_stop("User must verify very large run.");
    }
}

} // namespace feff::path
