// Read feff.inp and populate FeffInput structure.
// Converted from: src/RDINP/rdinp_l.f
// Original coded by S. Zabinski 1994, modified by A.L. Ankudinov March 2001.

#include "rdinp.hpp"
#include "iniall.hpp"
#include "setedg.hpp"
#include "ffsort.hpp"
#include "../common/string_utils.hpp"
#include "../common/itoken.hpp"
#include "../common/logging.hpp"
#include "../common/physics_utils.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <array>

namespace feff::rdinp {

// Helper: validate potential index range
static void phstop(int iph, const std::string& line) {
    if (iph < 0 || iph > nphx) {
        std::ostringstream ss;
        ss << " Unique potential index " << iph << " out of range."
           << " Must be between 0 and " << nphx << ".  Input line: " << line;
        feff::common::logger().wlog(ss.str());
        throw std::runtime_error("RDINP - PHSTOP");
    }
}

// Helper: print expert warning
static void warnex(const std::string& str) {
    auto& log = feff::common::logger();
    log.wlog(str);
    log.wlog(" Expert option, please read documentation carefully and check your results.");
}

int rdinp(FeffInput& inp, const std::string& input_file) {
    auto& log = feff::common::logger();

    // Constants
    constexpr int nssx = 16;
    constexpr int nbr = 30;
    constexpr double big = 1.0e5;

    // Initialize all input
    int nabs = 1;
    iniall(inp);

    // Local variables
    int iatom = 0;
    int ifolp = 0;
    int iovrlp = 0;
    int lxnat = 0;
    double folpx = 1.15;
    bool nogeom = false;
    inp.rclabs = big;
    double rmult = 1.0;
    double s02h = 1.0;
    inp.natt = 0;
    int ipr3 = 0;

    // Single scattering path data
    int nss = 0;
    std::array<int, nssx> indss{};
    std::array<int, nssx> iphss{};
    std::array<double, nssx> degss{};
    std::array<double, nssx> rss{};

    // Per-potential model atom indices (local)
    std::array<int, nphx + 1> iatph{};

    int icnt = 0; // EELS line counter

    // Open input file
    std::ifstream fin(input_file);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open input file: " + input_file);
    }

    // mode: 0=keyword, 1=atoms, 2=overlap, 3=potentials, 4=EELS
    int mode = 0;
    std::string line;

    while (std::getline(fin, line)) {
        // Strip carriage return (Windows line endings)
        if (!line.empty() && line.back() == '\r') line.pop_back();

        // Handle tabs and trim
        feff::common::replace_tabs(line);

        // Skip comments and blank lines
        if (feff::common::is_comment(line)) continue;

        // Split into words
        auto words = feff::common::split_words(line, 20);
        if (words.empty()) continue;

        int nwords = static_cast<int>(words.size());
        int itok = feff::common::itoken(words[0], input_file);

        // Label 210: process the card using current mode
    process_card:

        if (mode == 0) {
            // ---- Keyword mode ----

            if (itok == 1) {
                // ATOM(S) - following lines are atom positions
                mode = 1;
                iatom++;

            } else if (itok == 2) {
                // HOLE holecode s02
                if (!feff::common::parse_int(words[1], inp.ihole))
                    goto parse_error;
                if (nwords > 2)
                    feff::common::parse_double(words[2], s02h);
                mode = 0;

            } else if (itok == 3) {
                // OVERLAP iph
                int iph = 0;
                if (!feff::common::parse_int(words[1], iph))
                    goto parse_error;
                phstop(iph, line);
                warnex(" OVERLAP:");
                mode = 2;
                iovrlp++;
                // Store iph for reading overlap lines (we need to track which iph)
                // Use a small trick: store in iatph temporarily
                iatph[nphx] = iph; // temporary storage

            } else if (itok == 4) {
                // CONTROL
                if (nwords == 5) {
                    // feff7 input file
                    feff::common::parse_int(words[1], inp.mpot);
                    inp.mphase = inp.mpot;
                    inp.mfms = inp.mpot;
                    feff::common::parse_int(words[2], inp.mpath);
                    feff::common::parse_int(words[3], inp.mfeff);
                    feff::common::parse_int(words[4], inp.mchi);
                } else {
                    // feff8 input file
                    if (!feff::common::parse_int(words[1], inp.mpot)) goto parse_error;
                    if (!feff::common::parse_int(words[2], inp.mphase)) goto parse_error;
                    if (!feff::common::parse_int(words[3], inp.mfms)) goto parse_error;
                    if (!feff::common::parse_int(words[4], inp.mpath)) goto parse_error;
                    if (!feff::common::parse_int(words[5], inp.mfeff)) goto parse_error;
                    if (!feff::common::parse_int(words[6], inp.mchi)) goto parse_error;
                }
                mode = 0;

            } else if (itok == 5) {
                // EXCHANGE ixc vr0 vi0 (ixc0)
                inp.vr0 = 0.0;
                inp.vi0 = 0.0;
                if (nwords > 1) feff::common::parse_int(words[1], inp.ixc);
                if (nwords > 2) feff::common::parse_double(words[2], inp.vr0);
                if (nwords > 3) feff::common::parse_double(words[3], inp.vi0);
                if (nwords > 4) feff::common::parse_int(words[4], inp.ixc0);
                if (inp.ixc >= 3) warnex(" EXCHANGE >= 3:");
                mode = 0;

            } else if (itok == 6) {
                // ION iph xion(iph)
                int iph = 0;
                if (!feff::common::parse_int(words[1], iph)) goto parse_error;
                phstop(iph, line);
                if (!feff::common::parse_double(words[2], inp.xion[iph])) goto parse_error;
                warnex(" ION:");
                mode = 0;

            } else if (itok == 7) {
                // TITLE title...
                inp.ntitle++;
                if (inp.ntitle <= nheadx) {
                    // Extract everything after "TITLE "
                    std::string title_text;
                    if (line.size() > 5) {
                        title_text = line.substr(5);
                        title_text = feff::common::ltrim(title_text);
                    }
                    if (inp.ntitle <= static_cast<int>(inp.title.size())) {
                        inp.title[inp.ntitle - 1] = title_text;
                    }
                } else {
                    log.wlog(" Too many title lines, title ignored");
                    log.wlog(" " + line.substr(0, 71));
                }
                mode = 0;

            } else if (itok == 8) {
                // FOLP iph folp
                ifolp = 1;
                int iph = 0;
                if (!feff::common::parse_int(words[1], iph)) goto parse_error;
                phstop(iph, line);
                if (!feff::common::parse_double(words[2], inp.folp[iph])) goto parse_error;
                warnex(" FOLP:");
                mode = 0;

            } else if (itok == 9) {
                // RPATH rmax
                double rmax_d = 0.0;
                if (!feff::common::parse_double(words[1], rmax_d)) goto parse_error;
                inp.rmax = static_cast<float>(rmax_d);

            } else if (itok == 10) {
                // DEBYE temp debye-temp (idwopt)
                if (!feff::common::parse_double(words[1], inp.tk)) goto parse_error;
                if (!feff::common::parse_double(words[2], inp.thetad)) goto parse_error;
                inp.idwopt = 0;
                if (nwords > 3) {
                    feff::common::parse_int(words[3], inp.idwopt);
                    if (inp.idwopt > 2) {
                        std::ostringstream ss;
                        ss << " Option idwopt=" << inp.idwopt << " is not available.";
                        log.wlog(ss.str());
                        log.wlog("...setting idwopt=2 to use RM.");
                    }
                }
                mode = 0;

            } else if (itok == 11) {
                // RMULTIPLIER rmult
                if (!feff::common::parse_double(words[1], rmult)) goto parse_error;
                mode = 0;

            } else if (itok == 12) {
                // SS index ipot deg rss
                nss++;
                if (nss > nssx) {
                    std::ostringstream ss;
                    ss << " Too many ss paths requested, max is " << nssx;
                    log.wlog(ss.str());
                    throw std::runtime_error("RDINP: too many SS paths");
                }
                if (!feff::common::parse_int(words[1], indss[nss - 1])) goto parse_error;
                if (!feff::common::parse_int(words[2], iphss[nss - 1])) goto parse_error;
                if (!feff::common::parse_double(words[3], degss[nss - 1])) goto parse_error;
                if (!feff::common::parse_double(words[4], rss[nss - 1])) goto parse_error;
                mode = 0;

            } else if (itok == 13) {
                // PRINT ipr1 ipr2 ipr3 ipr4 ipr5 ipr6
                if (nwords == 5) {
                    // feff7 input
                    feff::common::parse_int(words[1], inp.ipr1);
                    inp.ipr2 = inp.ipr1;
                    ipr3 = inp.ipr1;
                    feff::common::parse_int(words[2], inp.ipr4);
                    feff::common::parse_int(words[3], inp.ipr5);
                    feff::common::parse_int(words[4], inp.ipr6);
                } else {
                    feff::common::parse_int(words[1], inp.ipr1);
                    feff::common::parse_int(words[2], inp.ipr2);
                    feff::common::parse_int(words[3], ipr3);
                    feff::common::parse_int(words[4], inp.ipr4);
                    feff::common::parse_int(words[5], inp.ipr5);
                    feff::common::parse_int(words[6], inp.ipr6);
                }
                mode = 0;

            } else if (itok == 14) {
                // POTENTIALS - following lines are potential definitions
                mode = 3;

            } else if (itok == 15) {
                // NLEG nlegmax
                if (!feff::common::parse_int(words[1], inp.nlegxx)) goto parse_error;
                mode = 0;

            } else if (itok == 16) {
                // CRITERIA critcw critpw
                if (!feff::common::parse_double(words[1], inp.critcw)) goto parse_error;
                double critpw_d = 0.0;
                if (!feff::common::parse_double(words[2], critpw_d)) goto parse_error;
                inp.critpw = static_cast<float>(critpw_d);
                mode = 0;

            } else if (itok == 17) {
                // NOGEOM
                nogeom = true;
                mode = 0;

            } else if (itok == 18) {
                // IORDER
                if (!feff::common::parse_int(words[1], inp.iorder)) goto parse_error;
                warnex(" IORDER:");
                mode = 0;

            } else if (itok == 19) {
                // PCRITERIA pcritk pcrith
                double pk = 0.0, ph = 0.0;
                if (!feff::common::parse_double(words[1], pk)) goto parse_error;
                if (!feff::common::parse_double(words[2], ph)) goto parse_error;
                inp.pcritk = static_cast<float>(pk);
                inp.pcrith = static_cast<float>(ph);
                mode = 0;

            } else if (itok == 20) {
                // SIG2 sig2g
                if (!feff::common::parse_double(words[1], inp.sig2g)) goto parse_error;
                mode = 0;

            } else if (itok == 21) {
                // XANES - not available
                log.wlog("XANES is not available in this version of feff8.5");
                log.wlog("Ignoring XANES");

            } else if (itok == 22) {
                // CORRECTIONS vrcorr vicorr
                if (!feff::common::parse_double(words[1], inp.vrcorr)) goto parse_error;
                if (!feff::common::parse_double(words[2], inp.vicorr)) goto parse_error;
                mode = 0;

            } else if (itok == 23) {
                // AFOLP (folpx)
                folpx = 1.15;
                if (nwords >= 2) feff::common::parse_double(words[1], folpx);
                mode = 0;

            } else if (itok == 24) {
                // EXAFS xkmax
                if (!feff::common::parse_double(words[1], inp.xkmax)) goto parse_error;
                mode = 0;

            } else if (itok == 25) {
                // POLARIZATION x y z
                inp.ipol = 1;
                if (!feff::common::parse_double(words[1], inp.evec[0])) goto parse_error;
                if (!feff::common::parse_double(words[2], inp.evec[1])) goto parse_error;
                if (!feff::common::parse_double(words[3], inp.evec[2])) goto parse_error;
                mode = 0;

            } else if (itok == 26) {
                // ELLIPTICITY elpty xivec(1:3)
                if (!feff::common::parse_double(words[1], inp.elpty)) goto parse_error;
                if (!feff::common::parse_double(words[2], inp.xivec[0])) goto parse_error;
                if (!feff::common::parse_double(words[3], inp.xivec[1])) goto parse_error;
                if (!feff::common::parse_double(words[4], inp.xivec[2])) goto parse_error;
                mode = 0;

            } else if (itok == 27) {
                // RGRID rgrd
                if (!feff::common::parse_double(words[1], inp.rgrd)) goto parse_error;
                warnex(" RGRID:");
                {
                    std::ostringstream ss;
                    ss << " RGRID, rgrd; " << std::scientific << std::setprecision(5) << inp.rgrd;
                    log.wlog(ss.str());
                }
                {
                    int i = 1 + static_cast<int>(12.5 / inp.rgrd);
                    if (i % 2 == 0) i++;
                    if (i > nrptx) {
                        std::ostringstream ss;
                        ss << " FATAL error in RGRID: increase nrptx to " << i;
                        log.wlog(ss.str());
                        throw std::runtime_error("RGRID: nrptx too small");
                    }
                }
                mode = 0;

            } else if (itok == 28) {
                // RPHASES
                warnex(" RPHASES:");
                log.wlog(" Real phase shifts only will be used.  FEFF results will be unreliable.");
                inp.lreal = 2;
                mode = 0;

            } else if (itok == 29) {
                // NSTAR
                inp.wnstar = true;
                mode = 0;

            } else if (itok == 30) {
                // NOHOLE
                if (inp.nohole < 0) {
                    inp.nohole = 0;
                    if (nwords >= 2) feff::common::parse_int(words[1], inp.nohole);
                    warnex(" NOHOLE:");
                }

            } else if (itok == 31) {
                // SIG3 alphat thetae
                if (!feff::common::parse_double(words[1], inp.alphat)) goto parse_error;
                if (nwords >= 3) feff::common::parse_double(words[2], inp.thetae);
                warnex(" SIG3:");
                {
                    std::ostringstream ss;
                    ss << " SIG3, alphat ; " << std::scientific << std::setprecision(5) << inp.alphat;
                    log.wlog(ss.str());
                }
                mode = 0;

            } else if (itok == 32) {
                // JUMPRM
                inp.jumprm = 1;

            } else if (itok == 33) {
                // MBCONV
                inp.mbconv = 1;

            } else if (itok == 34) {
                // SPIN ispin (spvec)
                if (!feff::common::parse_int(words[1], inp.ispin)) goto parse_error;
                if (inp.ispin != 0) inp.spvec[2] = 1.0;
                if (nwords > 2) feff::common::parse_double(words[2], inp.spvec[0]);
                if (nwords > 3) feff::common::parse_double(words[3], inp.spvec[1]);
                if (nwords > 4) feff::common::parse_double(words[4], inp.spvec[2]);

            } else if (itok == 35) {
                // EDGE L3
                setedg(words[1], inp.ihole);
                mode = 0;

            } else if (itok == 36) {
                // SCF rfms1 (lfms1 nscmt ca1 nmix ecv icoul)
                inp.nscmt = nbr;
                inp.ca1 = 0.2;
                double rfms1_d = 0.0;
                if (!feff::common::parse_double(words[1], rfms1_d)) goto parse_error;
                inp.rfms1 = static_cast<float>(rfms1_d);
                if (nwords > 2) feff::common::parse_int(words[2], inp.lfms1);
                if (nwords > 3) feff::common::parse_int(words[3], inp.nscmt);
                if (nwords > 4) feff::common::parse_double(words[4], inp.ca1);
                if (nwords > 5) feff::common::parse_int(words[5], inp.nmix);
                if (nwords > 6) feff::common::parse_double(words[6], inp.ecv);
                if (nwords > 7) feff::common::parse_int(words[7], inp.icoul);
                if (inp.nscmt <= 0 || inp.nscmt > nbr) inp.nscmt = nbr;
                if (inp.lfms1 > 0) inp.lfms1 = 1;
                if (inp.ca1 < 0.0) inp.ca1 = 0.0;
                if (inp.ca1 > 0.5) {
                    log.wlog(" Reduce convergence factors in SCF ");
                    throw std::runtime_error("Cannot run with specified ca1 in SCF card.");
                }
                if (inp.ecv >= 0.0) inp.ecv = -40.0;
                if (inp.nmix <= 0) inp.nmix = 1;
                if (inp.nmix > 30) inp.nmix = 30;

            } else if (itok == 37) {
                // FMS - not available
                log.wlog("FMS is not available in this version of FEFF8.5");
                log.wlog("Ignoring FMS");

            } else if (itok == 38) {
                // LDOS - not available
                log.wlog("LDOS is not available in this version of FEFF8.5");
                log.wlog("Ignoring LDOS");

            } else if (itok == 39) {
                // INTERSTITIAL inters (totvol)
                if (!feff::common::parse_int(words[1], inp.inters)) goto parse_error;
                if (nwords >= 3) feff::common::parse_double(words[2], inp.totvol);

            } else if (itok == 40) {
                // CFAVERAGE iphabs nabs rclabs
                if (!feff::common::parse_int(words[1], inp.iphabs)) goto parse_error;
                if (!feff::common::parse_int(words[2], nabs)) goto parse_error;
                double rclabs_d = 0.0;
                if (!feff::common::parse_double(words[3], rclabs_d)) goto parse_error;
                inp.rclabs = rclabs_d;
                if (inp.rclabs < 0.5) inp.rclabs = big;
                mode = 0;

            } else if (itok == 41) {
                // S02
                if (!feff::common::parse_double(words[1], inp.s02)) goto parse_error;
                mode = 0;

            } else if (itok == 42) {
                // XES - not available
                log.wlog("XES is not available in this version of FEFF8.5");
                log.wlog("Ignoring XES");

            } else if (itok == 43) {
                // DANES - not available
                log.wlog("DANES is not available in this version of FEFF8.5");
                log.wlog("Ignoring DANES");

            } else if (itok == 44) {
                // FPRIME - not available
                log.wlog("FPRIME is not available in this version of FEFF8.5");
                log.wlog("Ignoring FPRIME");

            } else if (itok == 45) {
                // RSIGMA
                warnex(" RSIGMA :");
                log.wlog(" Real self energy only will be used.  FEFF results will be unreliable.");
                if (inp.lreal < 1) inp.lreal = 1;
                mode = 0;

            } else if (itok == 46) {
                // XNCD/XMCD - not available
                log.wlog("XMCD is not available in this version of FEFF8.5");
                log.wlog("Ignoring XMCD");

            } else if (itok == 47) {
                // MULTIPOLES le2 (l2lp)
                if (!feff::common::parse_int(words[1], inp.le2)) goto parse_error;
                if (nwords > 2) feff::common::parse_int(words[2], inp.l2lp);
                mode = 0;

            } else if (itok == 48) {
                // UNFREEZEF
                inp.iunf = 1;
                mode = 0;

            } else if (itok == 49) {
                // TDLDA - not available
                log.wlog("TDLDA is not available in this version of FEFF8.5");
                log.wlog("Ignoring TDLDA");

            } else if (itok == 50) {
                // PMBSE - not available
                log.wlog("PMBSE is not available in this version of FEFF8.5");
                log.wlog("Ignoring PMBSE");

            } else if (itok == 51) {
                // PLASMON
                if (nwords > 1) {
                    feff::common::parse_int(words[1], inp.iPlsmn);
                } else {
                    inp.iPlsmn = 1;
                }

            } else if (itok == 52) {
                // S02CONV
                inp.mso2conv = 1;

            } else if (itok == 53) {
                // SELF
                inp.ipse = 1;

            } else if (itok == 54) {
                // SFSE k0
                inp.ipsk = 1;
                if (!feff::common::parse_double(words[1], inp.wsigk)) goto parse_error;

            } else if (itok == 55) {
                // RCONV cen cfname
                if (!feff::common::parse_double(words[1], inp.cen)) goto parse_error;
                inp.cfname = words[2].substr(0, 12);

            } else if (itok == 56) {
                // ELNES - not available
                log.wlog("ELNES is not available in this version of FEFF8.5");
                log.wlog("Ignoring ELNES");

            } else if (itok == 57) {
                // EXELFS - not available
                log.wlog("EXELFS is not available in this version of FEFF8.5");
                log.wlog("Ignoring EXELFS");

            } else if (itok == 58) {
                // MAGIC - not available
                log.wlog("MAGIC is not available in this version of FEFF8.5");
                log.wlog("Ignoring MAGIC");

            } else if (itok == 59) {
                // ABSOLUTE - not available
                log.wlog("ABSOLUTE is not available in this version of FEFF8.5");
                log.wlog("Ignoring ABSOLUTE");

            } else if (itok == 60) {
                // EGRID - not available
                log.wlog("EGRID is not available in this version of FEFF8.5");
                log.wlog("Ignoring EGRID");

            } else if (itok == -1) {
                // END
                break;

            } else {
                // Unrecognized keyword
                log.wlog(" " + line.substr(0, 70));
                log.wlog(" " + words[0]);
                {
                    std::ostringstream ss;
                    ss << " Token " << itok;
                    log.wlog(ss.str());
                }
                log.wlog(" Keyword unrecognized.");
                log.wlog(" See FEFF document -- some old features");
                log.wlog(" are no longer available.");
                throw std::runtime_error("RDINP-2: unrecognized keyword");
            }

        } else if (mode == 1) {
            // ---- Reading atom positions ----
            if (itok != 0) {
                // Done reading atoms, process current card as keyword
                mode = 0;
                goto process_card;
            }
            inp.natt++;
            if (inp.natt > nattx) {
                std::ostringstream ss;
                ss << "Too many atoms, max is " << nattx;
                log.wlog(ss.str());
                throw std::runtime_error("RDINP-3: too many atoms");
            }
            int idx = inp.natt - 1;
            if (!feff::common::parse_double(words[0], inp.ratx[idx][0])) goto parse_error;
            if (!feff::common::parse_double(words[1], inp.ratx[idx][1])) goto parse_error;
            if (!feff::common::parse_double(words[2], inp.ratx[idx][2])) goto parse_error;
            if (!feff::common::parse_int(words[3], inp.iphatx[idx])) goto parse_error;
            if (iatph[inp.iphatx[idx]] <= 0) iatph[inp.iphatx[idx]] = inp.natt;

        } else if (mode == 2) {
            // ---- Reading overlap instructions ----
            if (itok != 0) {
                mode = 0;
                goto process_card;
            }
            int iph = iatph[nphx]; // retrieve stored iph
            inp.novr[iph]++;
            int iovr_idx = inp.novr[iph] - 1;
            if (iovr_idx >= novrx) {
                std::ostringstream ss;
                ss << "Too many overlap shells, max is " << novrx;
                log.wlog(ss.str());
                throw std::runtime_error("RDINP-5: too many overlap shells");
            }
            if (!feff::common::parse_int(words[0], inp.iphovr[iph][iovr_idx])) goto parse_error;
            if (!feff::common::parse_int(words[1], inp.nnovr[iph][iovr_idx])) goto parse_error;
            if (!feff::common::parse_double(words[2], inp.rovr[iph][iovr_idx])) goto parse_error;

        } else if (mode == 3) {
            // ---- Reading unique potential definitions ----
            if (itok != 0) {
                mode = 0;
                goto process_card;
            }
            int iph = 0;
            if (!feff::common::parse_int(words[0], iph)) goto parse_error;
            if (iph < 0 || iph > nphx) {
                std::ostringstream ss;
                ss << "Unique potentials must be between 0 and " << nphx;
                log.wlog(ss.str());
                ss.str(""); ss.clear();
                ss << iph << " not allowed";
                log.wlog(ss.str());
                log.wlog(" " + line.substr(0, 71));
                throw std::runtime_error("RDINP: invalid potential index");
            }
            if (!feff::common::parse_int(words[1], inp.iz[iph])) goto parse_error;

            // Set default lmaxsc based on Z
            if (inp.iz[iph] < 6) {
                inp.lmaxsc[iph] = 1;
            } else if (inp.iz[iph] < 55) {
                inp.lmaxsc[iph] = 2;
            } else {
                inp.lmaxsc[iph] = 3;
            }

            // Optional potential label
            if (nwords >= 3) inp.potlbl[iph] = words[2].substr(0, 6);
            if (nwords >= 4) {
                int ltmp = 0;
                feff::common::parse_int(words[3], ltmp);
                if (ltmp >= 1 && ltmp <= lx) inp.lmaxsc[iph] = ltmp;
            }
            inp.lmaxph[iph] = 3;
            if (inp.iz[iph] < 6) inp.lmaxph[iph] = 2;
            if (nwords >= 5) {
                int ltmp = 0;
                feff::common::parse_int(words[4], ltmp);
                if (ltmp >= 1 && ltmp <= lx) inp.lmaxph[iph] = ltmp;
            }
            if (nwords >= 6) {
                feff::common::parse_double(words[5], inp.xnatph[iph]);
                lxnat = 1;
            }
            if (nwords >= 7) {
                feff::common::parse_double(words[6], inp.spinph[iph]);
            }

        } else if (mode == 4) {
            // ---- Reading EELS input (not available in this version) ----
            // Mode 4 would handle ELNES multi-line reading; skipped since
            // ELNES/EXELFS are disabled in feff8.5 EXAFS-only build
            mode = 0;

        } else {
            std::ostringstream ss;
            ss << "Mode unrecognized, mode " << mode;
            log.wlog(ss.str());
            throw std::runtime_error("RDINP-6: bad mode");
        }

        continue;

    parse_error:
        log.wlog(" Error reading input, bad line follows:");
        log.wlog(" " + line.substr(0, 71));
        throw std::runtime_error("RDINP fatal error: parse error");
    }

    fin.close();

    // ========================================================================
    // Post-processing: fix up defaults, error checks, etc.
    // ========================================================================

    // EELS checks (MAGIC, EELS, etc.) -- mostly no-ops for EXAFS-only build
    if (inp.magic == 1 && inp.ispec != 1) {
        log.wlog("To use MAGIC card you must have ELNES card.");
        log.wlog("Ignoring MAGIC card.");
        inp.magic = 0;
    }

    if (inp.ispec == 1 && inp.aver == 1 && inp.cross == 1) {
        log.wlog("WARNING: orientation averaged spectrum with cross-terms requested.");
        log.wlog("Averaging kills cross terms. Program ignores cross terms.");
    }

    // Set up EELS variables
    if (inp.ispec == 1) {
        if (inp.aver == 1) {
            inp.ipr1 = 1;
            inp.ipr1 = 10;
            inp.ipr1 = 10;
        } else {
            inp.ipr1 = 1;
            inp.ipr1 = 9;
            if (inp.cross == 1) {
                inp.ipr1 = 1;
            } else {
                inp.ipr1 = 4;
            }
        }
    }

    // Need smaller rgrid for nonlocal exchange
    if (inp.ixc0 < 0) inp.ixc0 = 0;
    if (inp.ixc % 10 >= 5 && inp.rgrd > 0.03) inp.rgrd = 0.03;
    if (inp.ixc0 % 10 >= 5 && inp.rgrd > 0.03) inp.rgrd = 0.03;

    // Must use linear polarization to use nstar
    if (inp.wnstar) {
        if (inp.ipol != 1) {
            log.wlog(" Must have linear polarization to use NSTAR.");
            log.wlog(" NSTAR will be turned off.");
            inp.wnstar = false;
        }
    }

    // Do not use ihole <= 0
    if (inp.ihole <= 0) {
        log.wlog(" Use NOHOLE to calculate without core hole.");
        log.wlog(" Only ihole greater than zero are allowed.");
        throw std::runtime_error("RDINP: ihole <= 0");
    }

    // Find how many unique potentials from POTENTIALS card
    inp.nph = 0;
    for (int iph = nphx; iph >= 0; --iph) {
        if (inp.iz[iph] > 0) {
            inp.nph = iph;
            break;
        }
    }

    // Cannot use OVERLAP and ATOMS together
    if (iatom > 0 && iovrlp > 0) {
        log.wlog(" Cannot use ATOMS and OVERLAP in the same feff.inp.");
        throw std::runtime_error("RDINP: ATOMS + OVERLAP conflict");
    }

    // Cannot use OVERLAP and CFAVERAGE together
    if (inp.novr[0] > 0) {
        inp.iphabs = 0;
        nabs = 1;
        inp.rclabs = big;
    }

    // Must have central atom (potential 0)
    if (inp.iz[0] <= 0) {
        if (inp.iphabs > 0) {
            inp.iz[0] = inp.iz[inp.iphabs];
            inp.potlbl[0] = inp.potlbl[inp.iphabs];
            inp.lmaxsc[0] = inp.lmaxsc[inp.iphabs];
            inp.lmaxph[0] = inp.lmaxph[inp.iphabs];
            inp.xion[0] = inp.xion[inp.iphabs];
        } else {
            log.wlog(" No absorbing atom (unique pot 0) was defined.");
            throw std::runtime_error("RDINP: no absorbing atom");
        }
    }

    // No gaps allowed in unique pots
    if (inp.iphabs > 0 && iatph[0] <= 0) iatph[0] = iatph[inp.iphabs];
    for (int iph = 0; iph <= inp.nph; ++iph) {
        if (iatph[iph] <= 0 && inp.novr[iph] <= 0) {
            std::ostringstream ss;
            ss << " No atoms or overlap cards for unique pot " << iph;
            log.wlog(ss.str());
            log.wlog(" Cannot calculate potentials, etc.");
            throw std::runtime_error("RDINP: missing atoms/overlap for potential");
        }
        // By default freeze f-electrons and reset lmaxsc=2
        if (inp.iunf == 0 && inp.lmaxsc[iph] > 2) inp.lmaxsc[iph] = 2;
    }

    // Count atoms per unique potential
    for (int iph = 0; iph <= inp.nph; ++iph) {
        if (lxnat == 0) {
            inp.xnatph[iph] = 0.0;
            for (int iat = 0; iat < inp.natt; ++iat) {
                if (inp.iphatx[iat] == iph) inp.xnatph[iph] += 1.0;
            }
            if (iph > 0 && iph == inp.iphabs)
                inp.xnatph[iph] -= 1.0;
        } else {
            if (inp.xnatph[iph] <= 0.01) {
                if (iph == 0) {
                    inp.xnatph[iph] = 0.01;
                } else {
                    std::ostringstream ss;
                    ss << " Inconsistency in POTENTIAL card is detected for unique pot " << iph;
                    log.wlog(ss.str());
                    log.wlog(" Results might be meaningless.");
                }
            }
        }
        if (inp.xnatph[iph] <= 0.0) inp.xnatph[iph] = 1.0;
    }

    // Normalize statistics if xnatph was given explicitly
    if (lxnat != 0) {
        for (int iph = 1; iph <= inp.nph; ++iph) {
            inp.xnatph[iph] = inp.xnatph[iph] / inp.xnatph[0];
        }
        inp.xnatph[0] = 1.0;
    }

    double xnat = 0.0;
    for (int iph = 0; iph <= inp.nph; ++iph) {
        xnat += inp.xnatph[iph];
    }

    // Find nearest and most distant atom distances
    double ratmin, ratmax;
    int icount = 0;
    if (inp.natt < 2) {
        ratmin = inp.rovr[0][0];
        ratmax = inp.rovr[0][inp.novr[0] - 1];
    } else {
        ratmax = 0.0;
        ratmin = 1.0e10;
        int iatabs = iatph[0];
        if (iatabs <= 0) iatabs = iatph[inp.iphabs];
        if (iatabs <= 0)
            throw std::runtime_error("RDINP fatal error: iatabs=NaN");

        for (int iat = 0; iat < inp.natt; ++iat) {
            if (inp.iphatx[iat] == inp.iphabs || inp.iphatx[iat] == 0)
                icount++;
            if (iat != (iatabs - 1)) {
                // Distance between atom iat and iatabs (convert to 0-based)
                double dx = inp.ratx[iat][0] - inp.ratx[iatabs - 1][0];
                double dy = inp.ratx[iat][1] - inp.ratx[iatabs - 1][1];
                double dz = inp.ratx[iat][2] - inp.ratx[iatabs - 1][2];
                double tmp = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (tmp > ratmax) ratmax = tmp;
                if (tmp < ratmin) ratmin = tmp;
            }
        }
        if (nabs <= 0) nabs = icount;
    }

    // Set total volume
    if (inp.totvol > 0.0)
        inp.totvol = inp.totvol * ratmin * ratmin * ratmin * xnat;

    // Set rfms if too small
    if (inp.rfms1 < ratmin) inp.rfms1 = -1.0f;
    if (inp.rfms2 < ratmin) inp.rfms2 = -1.0;
    if (inp.rfms2 < ratmin && inp.ispec < 2) inp.ispec = -inp.ispec;
    if (inp.rfms2 < ratmin && inp.ispec == 3) inp.ispec = -inp.ispec;

    // Set rmax if necessary
    if (inp.rmax <= 0.0f && nss <= 0 && inp.ispec <= 0) {
        inp.rmax = static_cast<float>(
            std::min(2.2 * ratmin, 1.01 * ratmax));
    }

    // Set core hole lifetime and s02
    feff::common::setgam(inp.iz[0], inp.ihole, inp.gamach);
    if (inp.s02 == 1.0) inp.s02 = s02h;

    // Apply rmult factor to distances (keep in Angstroms)
    inp.rmax = inp.rmax * static_cast<float>(rmult);
    inp.rfms1 = inp.rfms1 * static_cast<float>(rmult);
    inp.rfms2 = inp.rfms2 * rmult;
    inp.totvol = inp.totvol * rmult * rmult * rmult;

    for (int iat = 0; iat < inp.natt; ++iat) {
        for (int i = 0; i < 3; ++i) {
            inp.ratx[iat][i] *= rmult;
        }
    }
    for (int iph = 0; iph <= inp.nph; ++iph) {
        for (int iovr_idx = 0; iovr_idx < inp.novr[iph]; ++iovr_idx) {
            inp.rovr[iph][iovr_idx] *= rmult;
        }
    }
    for (int iss = 0; iss < nss; ++iss) {
        rss[iss] *= rmult;
    }

    // Clean up control flags
    if (inp.mpot != 0) inp.mpot = 1;
    if (inp.mphase != 0) inp.mphase = 1;
    if (inp.mfms != 0) inp.mfms = 1;
    if (inp.mpath != 0) inp.mpath = 1;
    if (inp.mfeff != 0) inp.mfeff = 1;
    if (inp.mchi != 0) inp.mchi = 1;
    if (nss <= 0) inp.ms = 1;
    if (ifolp != 0) inp.iafolp = -1;

    if (inp.natt <= 0) {
        // Overlap geometry only
        inp.mfms = 0;
        inp.mpath = 0;
        inp.ms = 0;
        inp.nscmt = 0;
        for (int iph = 0; iph <= inp.nph; ++iph) {
            if (inp.novr[iph] <= 0)
                throw std::runtime_error("Bad OVERLAP cards.");
        }
    }

    if (inp.iafolp >= 0) {
        for (int i = 0; i <= nphx; ++i) {
            inp.folp[i] = folpx;
        }
    }

    if (inp.ntitle <= 0) {
        inp.ntitle = 1;
        inp.title[0] = "Null title";
    }

    // Write SS paths.dat if needed
    if (nss > 0 && inp.mpath == 1) {
        std::ofstream fout("paths.dat");
        if (fout.is_open()) {
            for (int i = 0; i < inp.ntitle; ++i) {
                fout << " " << inp.title[i] << "\n";
            }
            fout << " Single scattering paths from ss lines cards in feff input\n";
            fout << " " << std::string(71, '-') << "\n";
            for (int iss = 0; iss < nss; ++iss) {
                if (inp.rmax <= 0.0f || rss[iss] <= inp.rmax) {
                    fout << std::setw(4) << indss[iss]
                         << std::setw(4) << 2
                         << std::fixed << std::setprecision(3) << std::setw(8) << degss[iss]
                         << "  index,nleg,degeneracy,r="
                         << std::setprecision(4) << std::setw(8) << rss[iss] << "\n";
                    fout << " single scattering\n";
                    fout << std::fixed << std::setprecision(6)
                         << std::setw(12) << rss[iss]
                         << std::setw(12) << 0.0
                         << std::setw(12) << 0.0
                         << std::setw(4) << iphss[iss]
                         << " '" << std::setw(6) << std::left
                         << inp.potlbl[iphss[iss]].substr(0, 6)
                         << std::right << "'\n";
                    fout << std::setw(12) << 0.0
                         << std::setw(12) << 0.0
                         << std::setw(12) << 0.0
                         << std::setw(4) << 0
                         << " '" << std::setw(6) << std::left
                         << inp.potlbl[0].substr(0, 6)
                         << std::right << "'  x,y,z,ipot\n";
                }
            }
            fout.close();
        }
    }

    // Log titles
    for (int i = 0; i < inp.ntitle; ++i) {
        log.wlog(" " + inp.title[i]);
    }

    // Build geometry
    if (nogeom) {
        if (inp.ipr4 < 2) inp.ipr4 = 2;
        if (nabs > 1)
            throw std::runtime_error("NOGEOM and CFAVERAGE are incompatible");
    } else {
        int iabs = 1;
        bool doptz = (inp.ispec != 1);
        ffsort(iabs, doptz, inp);
    }

    return nabs;
}

} // namespace feff::rdinp
