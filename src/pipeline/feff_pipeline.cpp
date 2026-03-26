// FeffPipeline — orchestrates the complete FEFF calculation.
// Converted from: src/feff6l/feff.f (main program)

#include "feff_pipeline.hpp"
#include "../par/parallel.hpp"
#include "../common/logging.hpp"
#include "../common/file_io.hpp"

// RDINP module
#include "../rdinp/iniall.hpp"
#include "../rdinp/rdinp.hpp"
#include "../rdinp/ffsort.hpp"

// JSON I/O
#include <feff/feff_input.hpp>
#include <feff/json_io.hpp>

// POT module
#include "../pot/reapot.hpp"
#include "../pot/pot.hpp"
#include "../pot/wrpot.hpp"

// XSPH module
#include "../xsph/rexsph.hpp"
#include "../xsph/xsph.hpp"

// PATH module
#include "../path/ffmod4.hpp"

// GENFMT module
#include "../genfmt/regenf.hpp"
#include "../genfmt/genfmt.hpp"

// FF2X module
#include "../ff2x/ff2x.hpp"

#include <iostream>
#include <stdexcept>
#include <cstring>

namespace feff::pipeline {

int run_feff(const PipelineConfig& config) {
    auto& log = feff::common::logger();

    try {
        // Initialize parallel environment
        feff::par::par_begin();

        log.open("feff.log");
        log.wlog("feff85exafs C++ version 8.5.0");
        log.wlog("");

        // Stage 1: Parse input
        log.wlog(" Reading input file...");
        int nabs = run_rdinp(config.input_file);
        if (nabs < 0) {
            log.wlog(" Error reading input file.");
            return 1;
        }

        // Stage 2: Potentials
        if (config.run_pot) {
            log.wlog(" Calculating potentials...");
            int rc = run_pot();
            if (rc != 0) {
                log.wlog(" Error in potential calculation.");
                return 2;
            }
        }

        // Stage 3: Phase shifts
        if (config.run_xsph) {
            log.wlog(" Calculating phase shifts...");
            int rc = run_xsph();
            if (rc != 0) {
                log.wlog(" Error in phase shift calculation.");
                return 3;
            }
        }

        // Stage 4: Path finding
        if (config.run_path) {
            log.wlog(" Finding scattering paths...");
            int rc = run_path();
            if (rc != 0) {
                log.wlog(" Error in path finding.");
                return 4;
            }
        }

        // Stage 5: F-matrix (GENFMT)
        if (config.run_genfmt) {
            log.wlog(" Calculating EXAFS parameters...");
            int rc = run_genfmt();
            if (rc != 0) {
                log.wlog(" Error in GENFMT calculation.");
                return 5;
            }
        }

        // Stage 6: Output spectra
        if (config.run_ff2x) {
            log.wlog(" Generating output spectra...");
            int rc = run_ff2x();
            if (rc != 0) {
                log.wlog(" Error generating output.");
                return 6;
            }
        }

        log.wlog("");
        log.wlog(" Feff done. Have a nice day.");

        feff::par::par_end();
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "FEFF error: " << e.what() << std::endl;
        log.wlog(std::string(" Fatal error: ") + e.what());
        return -1;
    }
}

// ============================================================================
// Stage 1: RDINP — Parse feff.inp, write JSON config files
// ============================================================================
int run_rdinp(const std::string& input_file) {
    auto& log = feff::common::logger();

    try {
        // 1. Create and initialize FeffInput with defaults
        feff::FeffInput inp;
        feff::rdinp::iniall(inp);

        // 2. Parse feff.inp — returns nabs (number of absorbers)
        log.wlog("  Parsing " + input_file + "...");
        int nabs = feff::rdinp::rdinp(inp, input_file);

        // 3. Write JSON interchange files for downstream modules
        log.wlog("  Writing JSON configuration files...");
        feff::json_io::write_pot_json(inp);
        feff::json_io::write_xsph_json(inp);
        feff::json_io::write_atoms_json(inp);
        feff::json_io::write_global_json(inp, nabs);
        feff::json_io::write_path_json(inp);
        feff::json_io::write_genfmt_json(inp);
        feff::json_io::write_ff2x_json(inp);

        // 4. Sort atoms and build geometry
        if (nabs == 1) {
            log.wlog("  Sorting atoms (single absorber)...");
            feff::rdinp::ffsort(1, true, inp);
        }

        // 5. Write geom.json from atom data for POT module
        //    Build the arrays that write_geom_json needs from inp
        {
            int nat = inp.natt;
            // Convert Vec3 array to double[][3]
            std::vector<std::array<double,3>> rat_arr(nat);
            std::vector<int> iphat_arr(nat);
            for (int i = 0; i < nat; ++i) {
                rat_arr[i] = {inp.ratx[i][0], inp.ratx[i][1], inp.ratx[i][2]};
                iphat_arr[i] = inp.iphatx[i];
            }
            // Build iatph: representative atom index per potential
            std::vector<int> iatph(inp.nph + 1, 0);
            for (int iph = 0; iph <= inp.nph; ++iph) {
                for (int i = 0; i < nat; ++i) {
                    if (iphat_arr[i] == iph) {
                        iatph[iph] = i;
                        break;
                    }
                }
            }
            feff::json_io::write_geom_json(inp, nat,
                reinterpret_cast<const double(*)[3]>(rat_arr.data()),
                iphat_arr.data(), iatph.data());
            log.wlog("  Wrote geom.json.");
        }

        log.wlog("  RDINP completed successfully.");
        return nabs;

    } catch (const std::exception& e) {
        log.wlog(std::string("  RDINP error: ") + e.what());
        return -1;
    }
}

// ============================================================================
// Stage 2: POT — Compute atomic potentials (SCF)
// ============================================================================
int run_pot() {
    auto& log = feff::common::logger();

    try {
        // --- Read input from JSON files via reapot ---
        int mpot = 0;
        double rgrd = 0.0;
        int ntitle = 0;
        std::string titles[nheadx];
        int ipr1 = 0, ispec = 0, nohole = 0, ihole = 0;
        double gamach = 0.0;
        int nph = 0;
        int iz[nphx + 1] = {};
        int lmaxsc[nphx + 1] = {};
        double xnatph[nphx + 1] = {};
        double xion[nphx + 1] = {};
        int iunf = 0, ixc = 0, jumprm = 0, iafolp = 0;
        double folp[nphx + 1] = {};
        int inters = 0;
        double totvol = 0.0;
        float rfms1 = 0.0f;
        int lfms1 = 0, nscmt = 0;
        double ca1 = 0.0;
        int nmix = 0;
        double ecv = 0.0;
        int icoul = 0;
        int novr[nphx + 1] = {};
        int iphovr_flat[(nphx + 1) * novrx] = {};
        int nnovr_flat[(nphx + 1) * novrx] = {};
        double rovr_flat[(nphx + 1) * novrx] = {};
        int nat = 0;
        double rat[3 * natx] = {};
        int iphat[natx] = {};
        int iatph[nphx + 1] = {};

        log.wlog("  Reading POT input from JSON...");
        feff::pot::reapot(mpot, rgrd, ntitle, titles,
                          ipr1, ispec, nohole, ihole, gamach,
                          nph, iz, lmaxsc, xnatph,
                          xion, iunf, ixc, jumprm, iafolp,
                          folp, inters, totvol,
                          rfms1, lfms1, nscmt, ca1, nmix,
                          ecv, icoul,
                          novr, iphovr_flat, nnovr_flat, rovr_flat,
                          nat, rat, iphat, iatph);

        // Check if potentials should be calculated
        if (mpot == 0) {
            log.wlog("  POT: mpot=0, skipping potential calculation.");
            return 0;
        }

        // --- Convert std::string titles to char[][80] for pot() ---
        char ctitle[nheadx][80];
        for (int i = 0; i < nheadx; ++i) {
            std::memset(ctitle[i], ' ', 80);
            auto len = std::min(titles[i].size(), static_cast<size_t>(79));
            std::memcpy(ctitle[i], titles[i].c_str(), len);
            ctitle[i][79] = '\0';
        }

        // --- Output arrays for pot() ---
        double rnrmav = 0.0, xmu = 0.0, vint = 0.0, rhoint = 0.0;
        double emu = 0.0, s02 = 0.0, erelax = 0.0, wp = 0.0;
        double rs = 0.0, xf = 0.0, qtotel = 0.0;
        int nhtmp = 0;

        int imt[nphx + 1] = {};
        double rmt[nphx + 1] = {};
        int inrm[nphx + 1] = {};
        double rnrm[nphx + 1] = {};
        double folpx[nphx + 1] = {};

        // Dirac spinor arrays
        double dgc0[251] = {};
        double dpc0[251] = {};

        constexpr int nph2 = nphx + 2;
        std::vector<double> dgc(251 * 30 * nph2, 0.0);
        std::vector<double> dpc_arr(251 * 30 * nph2, 0.0);
        std::vector<double> adgc(10 * 30 * nph2, 0.0);
        std::vector<double> adpc(10 * 30 * nph2, 0.0);

        // Density/potential arrays
        constexpr int nph1 = nphx + 1;
        std::vector<double> edens(251 * nph1, 0.0);
        std::vector<double> vclap(251 * nph1, 0.0);
        std::vector<double> vtot(251 * nph1, 0.0);
        std::vector<double> edenvl(251 * nph1, 0.0);
        std::vector<double> vvalgs(251 * nph1, 0.0);
        std::vector<double> dmag(251 * nph2, 0.0);
        std::vector<double> xnval(30 * nph2, 0.0);

        double eorb[30 * nph2] = {};
        int kappa[30 * nph2] = {};
        int iorb[8 * nph2] = {};
        double qnrm[nphx + 1] = {};
        double xnmues[(lx + 1) * nph1] = {};

        // --- Run the main pot calculation ---
        log.wlog("  Running POT kernel...");
        feff::pot::pot(true, rgrd, nohole,
                       inters, totvol, ecv, nscmt, nmix,
                       ntitle, ctitle,
                       nat, nph, ihole, iafolp, ixc,
                       iphat, rat, iatph, xnatph,
                       novr, iphovr_flat, nnovr_flat, rovr_flat,
                       folp, xion, iunf, iz, ipr1,
                       ispec, jumprm, lmaxsc, icoul,
                       ca1, rfms1, lfms1,
                       // output scalars
                       rnrmav, xmu, vint, rhoint,
                       emu, s02, erelax, wp,
                       rs, xf, qtotel,
                       // output arrays
                       imt, rmt, inrm, rnrm, folpx,
                       dgc0, dpc0,
                       dgc.data(), dpc_arr.data(), adgc.data(), adpc.data(),
                       edens.data(), vclap.data(), vtot.data(),
                       edenvl.data(), vvalgs.data(), dmag.data(), xnval.data(),
                       eorb, kappa, iorb,
                       qnrm, xnmues, nhtmp);

        // --- Write pot.pad ---
        log.wlog("  Writing pot.pad...");
        fprintf(stderr, "WRPOT_DEBUG: xmu=%.15e emu=%.15e erelax=%.15e s02=%.15e wp=%.15e\n",
                xmu, emu, erelax, wp, s02);
        feff::pot::wrpot(nph, ntitle, ctitle,
                         rnrmav, xmu, vint, rhoint,
                         emu, s02, erelax, wp,
                         ecv, rs, xf, qtotel,
                         imt, rmt, inrm, rnrm,
                         folp, folpx, xnatph,
                         dgc0, dpc0,
                         dgc.data(), dpc_arr.data(),
                         adgc.data(), adpc.data(),
                         edens.data(), vclap.data(), vtot.data(),
                         edenvl.data(), vvalgs.data(), dmag.data(),
                         xnval.data(),
                         eorb, kappa, iorb,
                         qnrm, xnmues,
                         nohole, ihole,
                         inters, totvol, iafolp,
                         xion, iunf, iz, jumprm);

        log.wlog("  POT completed successfully.");
        return 0;

    } catch (const std::exception& e) {
        log.wlog(std::string("  POT error: ") + e.what());
        return 1;
    }
}

// ============================================================================
// Stage 3: XSPH — Compute phase shifts and cross-sections
// ============================================================================
int run_xsph() {
    auto& log = feff::common::logger();

    try {
        // --- Read XSPH configuration from JSON ---
        feff::xsph::XsphConfig cfg;
        log.wlog("  Reading XSPH configuration...");
        feff::xsph::rexsph(cfg);

        // Check if phases should be calculated
        if (cfg.mphase == 0) {
            log.wlog("  XSPH: mphase=0, skipping phase shift calculation.");
            return 0;
        }

        // --- Read pot.pad for potential data ---
        log.wlog("  Reading pot.pad...");
        auto potdata = feff::common::read_pot("pot.pad");

        // --- Convert std::string titles to char[][80] for xsph() ---
        char ctitle[nheadx][80];
        for (int i = 0; i < nheadx; ++i) {
            std::memset(ctitle[i], ' ', 80);
            auto len = std::min(potdata.title[i].size(), static_cast<size_t>(79));
            std::memcpy(ctitle[i], potdata.title[i].c_str(), len);
            ctitle[i][79] = '\0';
        }

        // --- Allocate mutable copies of arrays that xsph may modify ---
        auto folpx = potdata.folpx;
        auto edens = potdata.edens;
        auto vclap = potdata.vclap;
        auto vtot = potdata.vtot;
        auto edenvl = potdata.edenvl;
        auto vvalgs = potdata.vvalgs;
        auto dmag = potdata.dmag;

        // --- Run XSPH calculation ---
        log.wlog("  Running XSPH kernel...");
        feff::xsph::xsph(
            true,                           // wrxsec: write xsect.json
            true,                           // verbose
            "phase.pad",                    // output file
            cfg.ipr2, cfg.ispec, cfg.vixan, cfg.xkstep, cfg.xkmax,
            cfg.gamach, cfg.rgrd,
            cfg.nph, cfg.lmaxph, cfg.potlbl,
            cfg.spinph, cfg.iatph, cfg.nat,
            &cfg.rat[0][0], cfg.iphat,
            cfg.ixc, cfg.vr0, cfg.vi0, cfg.ixc0, cfg.lreal,
            cfg.rfms2, cfg.lfms2, cfg.l2lp,
            cfg.ipol, cfg.ispin, cfg.le2, cfg.angks, cfg.ptz,
            cfg.iPl, cfg.izstd, cfg.ifxc, cfg.ipmbse, cfg.itdlda, cfg.nonlocal,
            potdata.ntitle, ctitle, potdata.rnrmav,
            potdata.xmu, potdata.vint, potdata.rhoint,
            potdata.emu, potdata.s02, potdata.erelax, potdata.wp, potdata.ecv,
            potdata.rs, potdata.xf, potdata.qtotel,
            potdata.imt.data(), potdata.rmt.data(),
            potdata.inrm.data(), potdata.rnrm.data(),
            potdata.folp.data(), folpx.data(),
            potdata.xnatph.data(),
            potdata.dgc0.data(), potdata.dpc0.data(),
            potdata.dgc.data(), potdata.dpc.data(),
            potdata.adgc.data(), potdata.adpc.data(),
            edens.data(), vclap.data(), vtot.data(),
            edenvl.data(), vvalgs.data(), dmag.data(),
            potdata.xnval.data(), potdata.iorb.data(),
            potdata.nohole, potdata.ihole,
            potdata.inters, potdata.totvol,
            potdata.iafolp,
            potdata.xion.data(), potdata.iunf, potdata.iz.data(), potdata.jumprm);

        log.wlog("  XSPH completed successfully.");
        return 0;

    } catch (const std::exception& e) {
        log.wlog(std::string("  XSPH error: ") + e.what());
        return 1;
    }
}

// ============================================================================
// Stage 4: PATH — Find scattering paths
// ============================================================================
int run_path() {
    auto& log = feff::common::logger();

    try {
        // The ffmod4 entry point handles the complete path-finding pipeline:
        //   repath() -> prcrit() -> paths() -> pathsd()
        log.wlog("  Running pathfinder (ffmod4)...");
        feff::path::ffmod4();

        log.wlog("  PATH completed successfully.");
        return 0;

    } catch (const std::exception& e) {
        log.wlog(std::string("  PATH error: ") + e.what());
        return 1;
    }
}

// ============================================================================
// Stage 5: GENFMT — Compute F-matrix scattering amplitudes
// ============================================================================
int run_genfmt() {
    auto& log = feff::common::logger();

    try {
        // --- Read genfmt configuration ---
        feff::genfmt::GenfmtConfig cfg;
        log.wlog("  Reading GENFMT configuration...");
        feff::genfmt::regenf(cfg);

        // Check if genfmt should run
        if (cfg.mfeff == 0) {
            log.wlog("  GENFMT: mfeff=0, skipping F-matrix calculation.");
            return 0;
        }

        // --- Run genfmt ---
        log.wlog("  Running GENFMT kernel...");
        feff::genfmt::genfmt(cfg.ipr5, cfg.critcw, cfg.iorder, cfg.wnstar,
                             cfg.ipol, cfg.ispin, cfg.le2, cfg.angks, cfg.elpty,
                             cfg.evec, cfg.xivec, cfg.ptz);

        log.wlog("  GENFMT completed successfully.");
        return 0;

    } catch (const std::exception& e) {
        log.wlog(std::string("  GENFMT error: ") + e.what());
        return 1;
    }
}

// ============================================================================
// Stage 6: FF2X — Generate output spectra (chi.dat, xmu.dat, etc.)
// ============================================================================
int run_ff2x() {
    auto& log = feff::common::logger();

    try {
        // The run_ff2x entry point in the ff2x module handles everything:
        //   read_ff2x_input() -> dispatch to ff2chi/ff2xmu/ff2afs based on ispec
        log.wlog("  Running FF2X module...");
        feff::ff2x::run_ff2x();

        log.wlog("  FF2X completed successfully.");
        return 0;

    } catch (const std::exception& e) {
        log.wlog(std::string("  FF2X error: ") + e.what());
        return 1;
    }
}

} // namespace feff::pipeline
