// Read XSPH module configuration from JSON files.
// Converted from src/XSPH/rexsph.f
//
// Reads from: geom.json, global.json, xsph.json
// Converts units to atomic units (Bohr, Hartree)

#include "rexsph.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/json_io.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <cstring>
#include <array>

namespace feff::xsph {

void rexsph(XsphConfig& config) {
    using json = nlohmann::json;

    // ---------------------------------------------------------------
    // 1. Read geometry from geom.json
    // ---------------------------------------------------------------
    {
        int ibounc[natx] = {};
        double rat3[natx][3] = {};
        int iphat[natx] = {};
        int iatph[nphx + 1] = {};

        feff::json_io::read_geom_json(config.nat, config.nph, iatph,
                                       rat3, iphat, ibounc);

        // Copy to config (transposed: Fortran rat(3,natx) → C++ rat[3][natx])
        for (int i = 0; i < config.nat; ++i) {
            config.rat[0][i] = rat3[i][0];
            config.rat[1][i] = rat3[i][1];
            config.rat[2][i] = rat3[i][2];
            config.iphat[i] = iphat[i];
        }
        for (int iph = 0; iph <= nphx; ++iph) {
            config.iatph[iph] = iatph[iph];
        }
    }

    // ---------------------------------------------------------------
    // 2. Read global params from global.json
    // ---------------------------------------------------------------
    {
        int nabs = 0, iphabs = 0;
        double rclabs = 0.0, elpty = 0.0;
        double evec[3] = {}, xivec[3] = {}, spvec[3] = {};
        std::array<std::array<FeffComplex, 3>, 3> ptz_arr{};

        feff::json_io::read_global_json(nabs, iphabs, rclabs,
                                         config.ipol, config.ispin, config.le2,
                                         elpty, config.angks,
                                         evec, xivec, spvec, ptz_arr);

        // Copy ptz to C-style array
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                config.ptz[i][j] = ptz_arr[i][j];
    }

    // ---------------------------------------------------------------
    // 3. Read xsph params from xsph.json
    // ---------------------------------------------------------------
    {
        std::ifstream ifs("xsph.json");
        if (!ifs.is_open()) {
            throw std::runtime_error("rexsph: cannot open xsph.json");
        }
        json j = json::parse(ifs);

        // Integer parameters
        config.mphase   = j.at("mphase").get<int>();
        config.ipr2     = j.at("ipr2").get<int>();
        config.ixc      = j.at("ixc").get<int>();
        config.ixc0     = j.at("ixc0").get<int>();
        config.ispec    = j.at("ispec").get<int>();
        config.lreal    = j.at("lreal").get<int>();
        config.lfms2    = j.at("lfms2").get<int>();
        config.nph      = j.at("nph").get<int>();
        config.l2lp     = j.at("l2lp").get<int>();
        config.iPl      = j.at("iPlsmn").get<int>();
        config.iGrid    = j.at("iGrid").get<int>();

        // Real parameters — note JSON keys match Fortran exactly
        // "vro" and "vio" (not "vr0"/"vi0") per the Fortran wrtjsn.f
        config.vr0      = j.at("vro").get<double>();
        config.vi0      = j.at("vio").get<double>();
        config.rgrd     = j.at("rgrd").get<double>();
        double rfms2_d  = j.at("rfms2").get<double>();
        config.rfms2    = static_cast<float>(rfms2_d);
        config.gamach   = j.at("gamach").get<double>();
        config.xkstep   = j.at("xkstep").get<double>();
        config.xkmax    = j.at("xkmax").get<double>();
        config.vixan    = j.at("vixan").get<double>();

        // Additional integer parameters
        config.izstd    = j.at("izstd").get<int>();
        config.ifxc     = j.at("ifxc").get<int>();
        config.ipmbse   = j.at("ipmbse").get<int>();
        config.itdlda   = j.at("itdlda").get<int>();
        config.nonlocal = j.at("nonlocal").get<int>();
        config.ibasis   = j.at("ibasis").get<int>();

        // Array parameters
        auto potlbl_arr = j.at("potlbl").get<std::vector<std::string>>();
        for (int i = 0; i < std::min(static_cast<int>(potlbl_arr.size()), nphx + 1); ++i) {
            std::strncpy(config.potlbl[i], potlbl_arr[i].c_str(), 6);
            config.potlbl[i][6] = '\0';
        }

        auto lmaxph_arr = j.at("lmaxph").get<std::vector<int>>();
        for (int i = 0; i < std::min(static_cast<int>(lmaxph_arr.size()), nphx + 1); ++i) {
            config.lmaxph[i] = lmaxph_arr[i];
        }

        auto spinph_arr = j.at("spinph").get<std::vector<double>>();
        for (int i = 0; i < std::min(static_cast<int>(spinph_arr.size()), nphx + 1); ++i) {
            config.spinph[i] = spinph_arr[i];
        }
    }

    // ---------------------------------------------------------------
    // 4. Convert to code units (Bohr and Hartrees)
    // ---------------------------------------------------------------
    config.rfms2  = config.rfms2 / static_cast<float>(feff::bohr);
    config.vr0    = config.vr0   / feff::hart;
    config.vi0    = config.vi0   / feff::hart;
    config.gamach = config.gamach / feff::hart;
    config.vixan  = config.vixan / feff::hart;
    config.xkstep = config.xkstep * feff::bohr;
    config.xkmax  = config.xkmax  * feff::bohr;

    for (int iat = 0; iat < config.nat; ++iat) {
        config.rat[0][iat] /= feff::bohr;
        config.rat[1][iat] /= feff::bohr;
        config.rat[2][iat] /= feff::bohr;
    }
}

} // namespace feff::xsph
