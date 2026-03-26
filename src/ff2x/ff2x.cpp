// Main FF2X module -- spectrum file naming, dispatch logic, and run_ff2x entry.
// Converted from: src/FF2X/ffmod6.f and shared utility code from ff2chi/ff2xmu/ff2afs.

#include "ff2x.hpp"
#include "reff2x.hpp"
#include "ff2chi.hpp"
#include "ff2xmu.hpp"
#include "ff2afs.hpp"

#include "../common/logging.hpp"

#include <stdexcept>
#include <string>

namespace feff::ff2x {

SpectrumFiles make_spectrum_files(int iip) {
    SpectrumFiles f;
    if (iip == 1) {
        f.chi_file  = "chi.dat";
        f.xmu_file  = "xmu.dat";
        f.pad_file  = "feff.pad";
        f.list_file = "list.dat";
    } else if (iip == 10) {
        f.chi_file  = "chi10.dat";
        f.xmu_file  = "xmu10.dat";
        f.pad_file  = "feff10.bin";
        f.list_file = "list10.dat";
    } else if (iip > 1 && iip < 10) {
        char c = static_cast<char>('0' + iip);
        f.chi_file  = std::string("chi0")  + c + ".dat";
        f.xmu_file  = std::string("xmu0")  + c + ".dat";
        f.pad_file  = std::string("feff0") + c + ".bin";
        f.list_file = std::string("list0") + c + ".dat";
    } else {
        throw std::runtime_error("Invalid spectrum index iip=" + std::to_string(iip));
    }
    return f;
}

bool is_cross_spectrum(int iip) {
    return !(iip == 1 || iip == 10 || iip == 5 || iip == 9);
}

void run_ff2x() {
    auto& log = common::logger();
    log.open("log6.dat");

    // Read input parameters
    FF2xParams p = read_ff2x_input();

    int iabs = 1;  // absorber index (configuration average loop handled externally)

    if (p.mchi == 1) {
        log.wlog(" Calculating chi...");

        if (p.ispec > 0 && p.ispec < 3) {
            // XANES: FMS+Paths method
            ff2xmu(p, iabs);
        } else if (p.ispec == 3 || p.ispec == 4) {
            // DANES/FPRIME: FMS+Paths method
            ff2afs(p, iabs);
        } else {
            // EXAFS: MS Paths expansion
            ff2chi(p, iabs);
        }

        log.wlog(" Done with module 6: DW + final sum over paths.");
    }

    log.close();
}

} // namespace feff::ff2x
