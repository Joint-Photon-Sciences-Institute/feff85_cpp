// Read FMS configuration from JSON files.
// Converted from: reafms.f
// Uses nlohmann/json to parse fms.json.
// Also reads geom.json and global.json via json_io module.

#include "reafms.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include <feff/constants.hpp>
#include <feff/json_io.hpp>

namespace feff::fms {

void reafms(FmsConfig& config,
            int& nat, int* iphat, double* rat,
            int& ipol, int& ispin, int& le2, double& angks,
            std::array<std::array<FeffComplex, 3>, 3>& ptz) {

    // Read geometry from geom.json
    int nph = 0;
    int iatph[nphx + 1] = {};
    int ibounc[natx] = {};

    feff::json_io::read_geom_json(nat, nph, iatph,
                                  reinterpret_cast<double(*)[3]>(rat),
                                  iphat, ibounc);

    // Read global parameters from global.json
    int nabs = 0, iphabs = 0;
    double rclabs = 0.0, elpty = 0.0;
    double evec[3] = {}, xivec[3] = {}, spvec[3] = {};

    feff::json_io::read_global_json(nabs, iphabs, rclabs,
                                    ipol, ispin, le2,
                                    elpty, angks,
                                    evec, xivec, spvec, ptz);

    // Read fms.json
    std::ifstream fin("fms.json");
    if (!fin.is_open()) {
        throw std::runtime_error("reafms: Failed to open fms.json");
    }

    nlohmann::json j;
    fin >> j;
    fin.close();

    auto get_or_throw = [&](const std::string& key) -> const nlohmann::json& {
        if (!j.contains(key)) {
            throw std::runtime_error("reafms: Missing key '" + key + "' in fms.json");
        }
        return j[key];
    };

    config.mfms   = get_or_throw("mfms").get<int>();
    config.idwopt = get_or_throw("idwopt").get<int>();
    config.minv   = get_or_throw("minv").get<int>();

    config.rfms2  = static_cast<float>(get_or_throw("rfms2").get<double>());
    config.rdirec = static_cast<float>(get_or_throw("rdirec").get<double>());
    config.toler1 = static_cast<float>(get_or_throw("toler1").get<double>());
    config.toler2 = static_cast<float>(get_or_throw("toler2").get<double>());

    config.tk     = get_or_throw("tk").get<double>();
    config.thetad = get_or_throw("thetad").get<double>();
    config.sig2g  = get_or_throw("sig2g").get<double>();

    auto lmaxph_arr = get_or_throw("lmaxph").get<std::vector<int>>();
    for (int iph = 0; iph <= nphx && iph < static_cast<int>(lmaxph_arr.size()); ++iph) {
        config.lmaxph[iph] = lmaxph_arr[iph];
    }

    // Convert from Angstrom to Bohr (matching Fortran reafms)
    config.rfms2  /= static_cast<float>(feff::bohr);
    config.rdirec /= static_cast<float>(feff::bohr);
    for (int iat = 0; iat < nat; ++iat) {
        for (int i = 0; i < 3; ++i) {
            rat[iat * 3 + i] /= feff::bohr;
        }
    }
}

} // namespace feff::fms
