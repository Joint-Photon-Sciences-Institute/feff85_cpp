// Input reading for FF2X module.
// Converted from: src/FF2X/reff2x.f
// Reads ff2x.json and global.json.

#include "reff2x.hpp"

#include <feff/constants.hpp>
#include <feff/json_io.hpp>

#include <nlohmann/json.hpp>
#include <fstream>
#include <array>
#include <complex>

namespace feff::ff2x {

FF2xParams read_ff2x_input() {
    FF2xParams p;

    // Read global.json (configuration average, polarization)
    int ipol = 0, ispin = 0, le2 = 0;
    double elpty = 0.0, angks = 0.0;
    double evec[3]{}, xivec[3]{}, spvec[3]{};
    std::array<std::array<std::complex<double>, 3>, 3> ptz{};

    json_io::read_global_json(p.nabs, p.iphabs, p.rclabs,
                              ipol, ispin, le2, elpty, angks,
                              evec, xivec, spvec, ptz);

    // Read ff2x.json
    // The json_io module should provide read_ff2x_json, but in the Fortran
    // code it was read directly with json_module. We use nlohmann/json.
    {
        nlohmann::json j;
        std::ifstream ifs("ff2x.json");
        if (!ifs.is_open()) {
            throw std::runtime_error("Failed to read ff2x.json");
        }
        ifs >> j;

        auto get_or_bail = [&](const char* key, auto& val) {
            if (!j.contains(key)) {
                json_io::bailout(key, "ff2x.json");
            }
            val = j[key].get<std::remove_reference_t<decltype(val)>>();
        };

        get_or_bail("mchi",    p.mchi);
        get_or_bail("ispec",   p.ispec);
        get_or_bail("idwopt",  p.idwopt);
        get_or_bail("ipr6",    p.ipr6);
        get_or_bail("mbconv",  p.mbconv);
        get_or_bail("absolu",  p.absolu);
        get_or_bail("vrcorr",  p.vrcorr);
        get_or_bail("vicorr",  p.vicorr);
        get_or_bail("s02",     p.s02);
        get_or_bail("critcw",  p.critcw);
        get_or_bail("tk",      p.tk);
        get_or_bail("thetad",  p.thetad);
        get_or_bail("alphat",  p.alphat);
        get_or_bail("thetae",  p.thetae);
        get_or_bail("sig2g",   p.sig2g);
    }

    // Default spectrum loop values
    p.elnes  = 0;
    p.ipstep = 1;
    p.ipmax  = 1;
    p.ipmin  = 1;

    // Convert energies from eV to code units (Hartree)
    p.vrcorr /= hart;
    p.vicorr /= hart;

    return p;
}

} // namespace feff::ff2x
