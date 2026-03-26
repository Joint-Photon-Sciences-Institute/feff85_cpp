// Main genfmt F-matrix generator.
// Converted from GENFMT/genfmt.f (~483 lines)
//
// Altered by Matt Newville (Jan 1999): format of feff.pad changed to packed-ascii.
// Altered by Alex Ankudinov (Feb 2000): disabled use of paths in LDOS.

#include "genfmt.hpp"
#include "genfmt_prep.hpp"
#include "rdpath.hpp"
#include "setlam.hpp"
#include "rot3i.hpp"
#include "sclmz.hpp"
#include "fmtrxi.hpp"
#include "mmtr.hpp"
#include "mmtrxi.hpp"
#include "import_calc.hpp"
#include "xstar.hpp"

#include "../common/logging.hpp"
#include "../common/pad_io.hpp"
#include "../common/physics_utils.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace feff::genfmt {

void genfmt(int ipr5, double critcw, int iorder, bool wnstar,
            int ipol, int ispin, int le2, double angks, double elpty,
            const double evec[3], const double xivec[3],
            const FeffComplex ptz[3][3]) {

    constexpr double eps = 1.0e-16;
    constexpr int mpadx = 8;

    // Data structures (replaces COMMON blocks)
    PhaseData pd;
    SpinPhaseData spd;
    NlmData nlm_data;
    RotationMatrixData rm;
    ClmzData clmz;
    FmatrixData fmat;
    LambdaData lam;
    BmatrixData bmatd;
    GenfmtPrepResult prep;

    // Local arrays
    FeffComplex rho[legtot];
    FeffComplex pmati[lamtot][lamtot][2];
    FeffComplex cchi[nex];
    FeffComplex rkk[nex][8];
    double xk[nex];
    double phff[nex], amff[nex];
    int lind[8];

    // ---- Initialization ----
    std::string phpad = "phase.pad";
    genfmt_prep(phpad, ispin, pd, spd, nlm_data, prep);

    int ne = prep.ne;
    int ne1 = prep.ne1;
    int npot = prep.npot;
    int lmaxp1 = prep.lmaxp1;
    int ilinit = prep.ilinit;
    int kinit = prep.kinit;

    int nsp = prep.nsp;
    int ll = prep.ll;

    // Copy xk from prep
    for (int ie = 0; ie < ne; ++ie) {
        xk[ie] = prep.xk[ie];
    }

    // ---- Open paths.dat and read header ----
    std::ifstream paths_in("paths.dat");
    common::check_file_open(paths_in, "paths.dat", "genfmt");

    // Read header lines (skip them) — mirrors Fortran rdhead.
    // Header ends at the "--------" separator line.
    std::string header_line;
    int ntext = 0;
    while (std::getline(paths_in, header_line)) {
        ntext++;
        // End-of-header marker: line starting with dashes
        if (header_line.find("--------") != std::string::npos) break;
    }

    // ---- Open list.dat ----
    std::ofstream list_out("list.dat");
    common::check_file_open(list_out, "list.dat", "genfmt");
    list_out << std::string(71, '-') << "\n";
    list_out << "  pathindex     sig2   amp ratio    "
             << "deg    nlegs  r effective\n";

    // ---- Open nstar.dat if necessary ----
    std::ofstream nstar_out;
    if (wnstar) {
        nstar_out.open("nstar.dat");
        nstar_out << " polarization" << std::fixed << std::setprecision(4)
                  << evec[0] << " " << evec[1] << " " << evec[2] << "\n";
        nstar_out << " npath     n*\n";
    }

    // ---- Importance threshold ----
    double crit0 = 0.0;
    if (ipr5 <= 0) crit0 = 2.0 * critcw / 3.0;

    {
        std::ostringstream ss;
        ss << "    Curved wave chi amplitude ratio" << std::fixed
           << std::setprecision(2) << std::setw(7) << critcw << "%";
        common::logger().wlog(ss.str());
    }
    if (ipr5 <= 0) {
        std::ostringstream ss;
        ss << "    Discard feff.dat for paths with cw ratio <"
           << std::fixed << std::setprecision(2) << std::setw(7) << crit0 << "%";
        common::logger().wlog(ss.str());
    }
    common::logger().wlog("    path  cw ratio     deg    nleg  reff");

    // ---- Open feff.pad ----
    std::ofstream feffpad("feff.pad");
    common::check_file_open(feffpad, "feff.pad", "genfmt");

    // Label line
    std::string vfeff_str = "Feff8L (EXAFS)";
    std::string vf85e_str = " 0.1";
    feffpad << "#_feff.pad v03: " << vfeff_str << vf85e_str << "\n";
    feffpad << "#_ " << npot << " " << ne << " " << mpadx << "\n";
    feffpad << "#& " << pd.ihole << " " << iorder << " " << ilinit
            << " " << std::scientific << std::setprecision(7)
            << pd.rnrmav << " " << pd.xmu << " " << pd.edge << "\n";

    // Pot labels and iz
    feffpad << "#@";
    for (int i = 0; i <= npot; ++i) {
        feffpad << " " << std::setw(6) << pd.potlbl[i];
    }
    for (int i = 0; i <= npot; ++i) {
        feffpad << " " << std::setw(3) << pd.iz[i];
    }
    feffpad << "\n";

    // Central atom phase shifts: ph(1:ne, ll, 0) as complex PAD
    {
        FeffComplex phc_tmp[nex];
        for (int ie = 0; ie < ne; ++ie) {
            phc_tmp[ie] = ph_at(pd, ie, ll, 0);
        }
        common::write_pad_complex(feffpad, mpadx, phc_tmp, ne);
    }
    // ck(1:ne) as complex PAD
    common::write_pad_complex(feffpad, mpadx, prep.ck, ne);
    // xkr(1:ne) as real PAD
    common::write_pad_double(feffpad, mpadx, prep.xkr, ne);

    // ---- Main loop over paths ----
    bool done = false;
    while (!done) {
        rdpath(paths_in, done, ipol, pd.potlbl, pd.rat, pd.ri, pd.beta, pd.eta,
               pd.deg, pd.ipot, pd.nsc, pd.nleg, npot, pd.ipath);

        if (done) break;

        int icalc = iorder;
        prep.npath++;
        prep.ntotal++;

        int nleg = pd.nleg;
        int nsc = pd.nsc;
        double deg = pd.deg;

        // N* star calculation
        if (wnstar) {
            double eps1[3], eps2[3], vec1[3], vec2[3];
            for (int ic = 0; ic < 3; ++ic) {
                vec1[ic] = pd.rat[ic][1] - pd.rat[ic][0];
                vec2[ic] = pd.rat[ic][nleg - 1] - pd.rat[ic][0];
                eps1[ic] = evec[ic];
            }
            double eps2_arr[3] = {};
            if (elpty != 0.0) {
                eps2_arr[0] = xivec[1] * evec[2] - xivec[2] * evec[1];
                eps2_arr[1] = xivec[2] * evec[0] - xivec[0] * evec[2];
                eps2_arr[2] = xivec[0] * evec[1] - xivec[1] * evec[0];
            }
            int ndeg = static_cast<int>(std::round(deg));
            double xxstar_val = xstar(eps1, eps2_arr, vec1, vec2, ndeg, elpty, ilinit);
            nstar_out << std::setw(6) << prep.npath << std::fixed
                      << std::setprecision(3) << std::setw(10) << xxstar_val << "\n";
        }

        // Compute reff
        double reff = 0.0;
        for (int i = 1; i <= nleg; ++i) {
            reff += pd.ri[i];
        }
        reff /= 2.0;

        // Set lambda for low k
        setlam(icalc, 1, pd.beta, nsc, nleg, ilinit, lam);

        // Calculate and store rotation matrix elements
        // Fortran rot3i uses 1-based ileg for both beta access and storage.
        // C++ fmtrxi uses 0-based ilegp. We pass isc (1-based) for beta
        // access, and store at isc-1 (0-based) for fmtrxi compatibility.
        for (int isc = 1; isc <= nleg; ++isc) {
            rot3i(lmaxp1, lam.mmaxp1, isc, pd.beta, rm, isc - 1);
        }
        if (ipol > 0) {
            rot3i(ilinit + 1, ilinit + 1, nleg + 1, pd.beta, rm, nleg);
        }

        // Zero cchi
        for (int ie = 0; ie < ne; ++ie) {
            cchi[ie] = FeffComplex(0.0, 0.0);
        }

        // ---- Spin loop ----
        for (int is = 0; is < nsp; ++is) {
            int ispin_eff = (nsp == 1) ? ispin : is + 1;
            mmtr(bmatd, ipol, ispin_eff, le2, angks, ptz, lind,
                 rm, pd.eta, nsc, nleg, kinit, ilinit);

            // Copy spin-specific data
            for (int ie = 0; ie < ne; ++ie) {
                pd.eref[ie] = spd.eref2[ie][is];
            }
            for (int iph = 0; iph <= npot; ++iph) {
                for (int ie = 0; ie < ne; ++ie) {
                    for (int il = -pd.lmax[ie][iph]; il <= pd.lmax[ie][iph]; ++il) {
                        ph_at(pd, ie, il, iph) = ph4_at(spd, ie, il, is, iph);
                    }
                }
            }
            for (int ie = 0; ie < ne; ++ie) {
                for (int kdif = 0; kdif < 8; ++kdif) {
                    rkk[ie][kdif] = spd.rkk2[ie][kdif][is];
                }
            }
            FeffComplex ck_local[nex];
            for (int ie = 0; ie < ne; ++ie) {
                ck_local[ie] = std::sqrt(2.0 * (pd.em[ie] - pd.eref[ie]));
            }

            // ---- Energy loop ----
            for (int ie = 0; ie < ne; ++ie) {
                // Complex rho
                for (int ileg = 1; ileg <= nleg; ++ileg) {
                    rho[ileg - 1] = ck_local[ie] * pd.ri[ileg];
                }

                if (std::abs(ck_local[ie]) <= eps) {
                    // ck is zero, xafs undefined
                    std::ostringstream ss;
                    ss << " genfmt: ck=0.  ie, ck(ie) " << ie
                       << " " << ck_local[ie].real() << " " << ck_local[ie].imag();
                    common::logger().wlog(ss.str());
                    continue;
                }

                // Zero clmi arrays
                for (int il = 0; il < legtot; ++il)
                    for (int im = 0; im < mtot + ntot + 1; ++im)
                        for (int il2 = 0; il2 < ltot + 1; ++il2)
                            clmz.clmi[il2][im][il] = FeffComplex(0.0, 0.0);

                int mnmxp1 = lam.mmaxp1 + lam.nmax;

                // Compute spherical wave factors for each leg
                // Note: sclmz uses 0-based ileg index
                for (int ileg = 1; ileg <= nleg; ++ileg) {
                    int isc0 = ileg - 1;
                    if (isc0 == 0) isc0 = nleg;
                    int isc1 = ileg;
                    int lxp1 = std::max(pd.lmax[ie][pd.ipot[isc0]] + 1,
                                        pd.lmax[ie][pd.ipot[isc1]] + 1);
                    int mnp1 = std::min(lxp1, mnmxp1);
                    sclmz(rho, lxp1, mnp1, ileg - 1, clmz);
                }

                // Calculate scattering matrices
                // First matrix (Fortran: ie=ie, ileg=2, ilegp=1)
                // C++ 0-based: ileg=1, ilegp=0
                fmtrxi(lam.lamx, lam.laml0x, ie, 1, 0,
                       clmz, lam, nlm_data, rm, pd, fmat);

                // Last matrix if nleg > 2
                if (nleg > 2) {
                    fmtrxi(lam.laml0x, lam.lamx, ie, nleg - 1, nleg - 2,
                           clmz, lam, nlm_data, rm, pd, fmat);
                }

                // Intermediate scattering matrices
                for (int ilegp = 2; ilegp <= nsc - 1; ++ilegp) {
                    int ileg_idx = ilegp;  // 0-based
                    fmtrxi(lam.lamx, lam.lamx, ie, ileg_idx, ilegp - 1,
                           clmz, lam, nlm_data, rm, pd, fmat);
                }

                // ---- Matrix multiplication (trace of product) ----
                // f(2,1) -> pmat(0)
                int indp = 0;
                for (int lmp = 0; lmp < lam.laml0x; ++lmp) {
                    for (int lm = 0; lm < lam.lamx; ++lm) {
                        pmati[lm][lmp][indp] = fmat.fmati[lm][lmp][0];
                    }
                }

                // f(N,N-1) * ... * f(3,2) * [f(2,1)]
                for (int isc = 2; isc <= nleg - 1; ++isc) {
                    indp = 1 - (isc % 2);  // Fortran: 2 - mod(isc,2)
                    int indp0 = 1 - indp;   // Fortran: 1 + mod(indp,2)

                    for (int lmp = 0; lmp < lam.laml0x; ++lmp) {
                        for (int lm = 0; lm < lam.lamx; ++lm) {
                            FeffComplex pllp(0.0, 0.0);
                            for (int lmi = 0; lmi < lam.lamx; ++lmi) {
                                pllp += fmat.fmati[lm][lmi][isc - 1] *
                                        pmati[lmi][lmp][indp0];
                            }
                            pmati[lm][lmp][indp] = pllp;
                        }
                    }
                }

                // srho = sum rho(i), prho = prod rho(i)
                FeffComplex srho(0.0, 0.0);
                FeffComplex prho(1.0, 0.0);
                for (int ileg = 0; ileg < nleg; ++ileg) {
                    srho += rho[ileg];
                    prho *= rho[ileg];
                }

                // Termination matrix
                mmtrxi(rkk, lam.laml0x, bmatd, ie, 0, nleg - 1, lind,
                       clmz, lam, nlm_data, pd.eta, fmat);

                // Final trace
                FeffComplex ptrac(0.0, 0.0);
                for (int lm = 0; lm < lam.laml0x; ++lm) {
                    for (int lmp = 0; lmp < lam.laml0x; ++lmp) {
                        ptrac += fmat.fmati[lm][lmp][nleg - 1] *
                                 pmati[lmp][lm][indp];
                    }
                }

                // Complex chi (without 2kr term)
                FeffComplex cfac = std::exp(coni * (srho - 2.0 * xk[ie] * reff)) / prho;
                if (nsp == 2 && is == 0) cfac = -cfac;
                cchi[ie] += ptrac * cfac;

            }
            // end energy loop
        }
        // end spin loop

        // Compute importance factor
        double crit;
        import_calc(ne1, nsp, prep.ik0, reff, deg, prep.ckmag,
                    pd.em, spd.eref2, cchi, prep.xportx, crit);

        // Write path data to feff.pad if important enough
        if (ipr5 >= 1 || crit >= crit0) {
            // Write path info to feff.pad
            feffpad << "##" << std::setw(6) << pd.ipath
                    << " " << std::setw(3) << nleg
                    << " " << std::fixed << std::setprecision(3) << std::setw(7) << deg
                    << " " << std::setprecision(7) << std::setw(11) << reff * bohr
                    << " " << std::scientific << std::setprecision(4) << std::setw(15) << crit;
            for (int i = 1; i <= nleg; ++i) {
                feffpad << " " << std::setw(2) << pd.ipot[i];
            }
            feffpad << "\n";

            // Write per-path geometry arrays in PAD format
            // rat(3*nleg): positions flattened as x1,y1,z1, x2,y2,z2, ...
            {
                std::vector<double> rat_flat(3 * nleg);
                for (int j = 0; j < nleg; ++j) {
                    rat_flat[3 * j]     = pd.rat[0][j + 1];
                    rat_flat[3 * j + 1] = pd.rat[1][j + 1];
                    rat_flat[3 * j + 2] = pd.rat[2][j + 1];
                }
                common::write_pad_double(feffpad, mpadx, rat_flat.data(), 3 * nleg);
            }
            // beta(nleg): Fortran writes beta(1:nleg)
            common::write_pad_double(feffpad, mpadx, &pd.beta[1], nleg);
            // eta(nleg): Fortran writes eta(0:nleg-1)
            common::write_pad_double(feffpad, mpadx, &pd.eta[0], nleg);
            // ri(nleg): Fortran writes ri(1:nleg)
            common::write_pad_double(feffpad, mpadx, &pd.ri[1], nleg);

            // Compute phase and amplitude arrays
            double phffo = 0.0;
            for (int ie = 0; ie < ne; ++ie) {
                phff[ie] = 0.0;
                if (std::abs(cchi[ie]) >= eps) {
                    phff[ie] = std::atan2(cchi[ie].imag(), cchi[ie].real());
                }
                if (ie > 0) {
                    common::pijump(phff[ie], phffo);
                }
                phffo = phff[ie];
                amff[ie] = std::abs(cchi[ie]);
            }

            // Write amplitude and phase arrays in PAD format
            common::write_pad_double(feffpad, mpadx, amff, ne);
            common::write_pad_double(feffpad, mpadx, phff, ne);

            // Write to list.dat
            list_out << std::setw(8) << pd.ipath
                     << std::fixed << std::setprecision(5) << std::setw(12) << 0.0
                     << std::setprecision(3) << std::setw(10) << crit
                     << std::setw(10) << deg
                     << std::setw(6) << nleg
                     << std::setprecision(4) << std::setw(9) << reff * bohr << "\n";

            // Tell user
            {
                std::ostringstream ss;
                ss << "   " << std::setw(4) << pd.ipath
                   << std::fixed << std::setprecision(3)
                   << std::setw(10) << crit << std::setw(10) << deg
                   << std::setw(6) << nleg
                   << std::setprecision(4) << std::setw(9) << reff * bohr;
                common::logger().wlog(ss.str());
            }
            prep.nused++;
        } else {
            std::ostringstream ss;
            ss << "   " << std::setw(4) << pd.ipath
               << std::fixed << std::setprecision(3)
               << std::setw(10) << crit << std::setw(10) << deg
               << std::setw(6) << nleg
               << std::setprecision(4) << std::setw(9) << reff * bohr
               << " neglected";
            common::logger().wlog(ss.str());
        }
    }
    // end path loop

    // Close files
    paths_in.close();
    list_out.close();
    feffpad.close();
    if (wnstar) nstar_out.close();

    {
        std::ostringstream ss;
        ss << " " << prep.nused << " paths kept, " << prep.ntotal << " examined.";
        common::logger().wlog(ss.str());
    }
}

} // namespace feff::genfmt
