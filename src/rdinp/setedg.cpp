// Map edge name to core-hole index.
// Converted from: src/RDINP/setedg.f

#include "setedg.hpp"
#include <stdexcept>
#include <array>

namespace feff::rdinp {

void setedg(const std::string& edge_name, int& ihole) {
    // Edge labels: spectroscopic notation
    static const std::array<const char*, 30> edglbl = {
        "NO", "K",  "L1", "L2", "L3",
        "M1", "M2", "M3", "M4", "M5",
        "N1", "N2", "N3", "N4", "N5", "N6", "N7",
        "O1", "O2", "O3", "O4", "O5", "O6", "O7",
        "P1", "P2", "P3", "P4", "P5", "R1"
    };

    // Edge labels: numeric notation ("0" through "29")
    static const std::array<const char*, 30> edglbp = {
        "0",  "1",  "2",  "3",  "4",
        "5",  "6",  "7",  "8",  "9",
        "10", "11", "12", "13", "14", "15", "16",
        "17", "18", "19", "20", "21", "22", "23",
        "24", "25", "26", "27", "28", "29"
    };

    ihole = -1;
    for (int i = 0; i < 30; ++i) {
        if (edge_name == edglbl[i] || edge_name == edglbp[i]) {
            ihole = i;
        }
    }

    if (ihole < 0) {
        throw std::runtime_error("unknown EDGE: " + edge_name);
    }
}

} // namespace feff::rdinp
