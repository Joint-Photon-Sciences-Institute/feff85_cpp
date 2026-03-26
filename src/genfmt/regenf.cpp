// Read genfmt configuration from genfmt.json and global.json.
// Converted from GENFMT/regenf.f
//
// Reads genfmt.json for module-specific settings (mfeff, ipr5, iorder, etc.)
// and global.json for polarization parameters (ipol, evec, ptz, etc.).

#include "regenf.hpp"
#include <feff/json_io.hpp>

#include <nlohmann/json.hpp>
#include <array>
#include <complex>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;

namespace feff::genfmt {

static json read_json(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("regenf: cannot open " + filename);
    }
    return json::parse(ifs);
}

void regenf(GenfmtConfig& config) {
    // ---- Read genfmt.json ----
    {
        json j = read_json("genfmt.json");
        config.mfeff  = j.value("mfeff",  1);
        config.ipr5   = j.value("ipr5",   0);
        config.iorder = j.value("iorder", 2);
        config.wnstar = j.value("wnstar", false);
        config.critcw = j.value("critcw", 4.0);
    }

    // ---- Read global.json for polarization parameters ----
    {
        int nabs = 1, iphabs = 0;
        double rclabs = 0.0;
        double spvec[3] = {};
        std::array<std::array<std::complex<double>, 3>, 3> ptz_std{};

        json_io::read_global_json(nabs, iphabs, rclabs,
                                  config.ipol, config.ispin, config.le2,
                                  config.elpty, config.angks,
                                  config.evec, config.xivec, spvec, ptz_std);

        // Convert ptz from std::array<std::array<complex,3>,3> to FeffComplex[3][3]
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                config.ptz[i][j] = FeffComplex(ptz_std[i][j].real(),
                                               ptz_std[i][j].imag());
            }
        }
    }
}

} // namespace feff::genfmt
