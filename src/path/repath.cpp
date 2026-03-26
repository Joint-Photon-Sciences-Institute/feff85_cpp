// Read input for the PATH module.
// Converted from: src/PATH/repath.f
// Uses JSON I/O instead of Fortran formatted reads.

#include "repath.hpp"
#include <feff/json_io.hpp>
#include <feff/dimensions.hpp>
#include <nlohmann/json.hpp>
#include <complex>
#include <array>
#include <fstream>
#include <stdexcept>

namespace feff::path {

void repath(int& ms, int& mpath, int& ipr4,
            float& pcritk, float& pcrith, int& nncrit,
            float& rmax, int& nlegxx, float& rfms2, float& critpw,
            int& nat, double rat[][3], int iphat[], int ibounc[],
            int& ipol, int& ispin, double evec[3], double xivec[3]) {

    // Read geometry
    int nph = 0;
    int iatph[nphx + 1] = {};
    feff::json_io::read_geom_json(nat, nph, iatph, rat, iphat, ibounc);

    // Read global parameters
    int nabs, iphabs, le2;
    double rclabs, elpty, angks;
    double spvec[3];
    std::array<std::array<std::complex<double>, 3>, 3> ptz;
    feff::json_io::read_global_json(nabs, iphabs, rclabs,
                                     ipol, ispin, le2,
                                     elpty, angks,
                                     evec, xivec, spvec, ptz);

    // Read path.json
    std::ifstream fin("path.json");
    if (!fin) {
        throw std::runtime_error("Failed to read path.json");
    }
    nlohmann::json j;
    fin >> j;

    mpath  = j.at("mpath").get<int>();
    ms     = j.at("ms").get<int>();
    nncrit = j.at("nncrit").get<int>();
    nlegxx = j.at("nlegxx").get<int>();
    ipr4   = j.at("ipr4").get<int>();
    critpw = static_cast<float>(j.at("critpw").get<double>());
    pcritk = static_cast<float>(j.at("pcritk").get<double>());
    pcrith = static_cast<float>(j.at("pcrith").get<double>());
    rmax   = static_cast<float>(j.at("rmax").get<double>());
    rfms2  = static_cast<float>(j.at("rfms2").get<double>());
}

} // namespace feff::path
