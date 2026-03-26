#pragma once
// Periodic table data: atomic weights and symbols for Z=1..103.
// Converted from: src/COMMON/pertab.f
// Data source: Ashcroft & Mermin (inside front cover)

#include <array>
#include <stdexcept>
#include <string>

namespace feff::common {

// Number of elements in the table (H=1 through Lr=103)
inline constexpr int num_elements = 103;

// Atomic weights indexed by Z (1-based; index 0 is unused/zero)
inline constexpr std::array<double, num_elements + 1> atomic_weights = {{
    0.0,       // placeholder for index 0
    1.0079,    4.0026,    6.941,     9.0122,    10.81,     12.01,     // H-C
    14.007,    15.999,    18.998,    20.18,     22.9898,   24.305,    // N-Mg
    26.982,    28.086,    30.974,    32.064,    35.453,    39.948,    // Al-Ar
    39.09,     40.08,     44.956,    47.90,     50.942,    52.00,     // K-Cr
    54.938,    55.85,     58.93,     58.71,     63.55,     65.38,     // Mn-Zn
    69.72,     72.59,     74.922,    78.96,     79.91,     83.80,     // Ga-Kr
    85.47,     87.62,     88.91,     91.22,     92.91,     95.94,     // Rb-Mo
    98.91,     101.07,    102.90,    106.40,    107.87,    112.40,    // Tc-Cd
    114.82,    118.69,    121.75,    127.60,    126.90,    131.30,    // In-Xe
    132.91,    137.34,    138.91,    140.12,    140.91,    144.24,    // Cs-Nd
    145.0,     150.35,    151.96,    157.25,    158.92,    162.50,    // Pm-Dy
    164.93,    167.26,    168.93,    173.04,    174.97,    178.49,    // Ho-Hf
    180.95,    183.85,    186.2,     190.20,    192.22,    195.09,    // Ta-Pt
    196.97,    200.59,    204.37,    207.19,    208.98,    210.0,     // Au-Po
    210.0,     222.0,     223.0,     226.0,     227.0,     232.04,    // At-Th
    231.0,     238.03,    237.05,    244.0,     243.0,     247.0,     // Pa-Cm
    247.0,     251.0,     254.0,     257.0,     256.0,     254.0,     // Bk-No
    257.0                                                              // Lr
}};

// Atomic symbols indexed by Z (1-based; index 0 is "?")
inline constexpr std::array<const char*, num_elements + 1> atomic_symbols = {{
    "?",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr"
}};

/// Get atomic weight for element with atomic number iz (1-103).
inline double atomic_weight(int iz) {
    if (iz < 1 || iz > num_elements) {
        throw std::out_of_range("atomic_weight: iz=" + std::to_string(iz) + " out of range [1,103]");
    }
    return atomic_weights[iz];
}

/// Get atomic symbol for element with atomic number iz (1-103).
inline const char* atomic_symbol(int iz) {
    if (iz < 1 || iz > num_elements) {
        throw std::out_of_range("atomic_symbol: iz=" + std::to_string(iz) + " out of range [1,103]");
    }
    return atomic_symbols[iz];
}

/// Look up atomic number from symbol string (case-insensitive). Returns 0 if not found.
inline int atomic_number(const std::string& sym) {
    for (int iz = 1; iz <= num_elements; ++iz) {
        std::string s = atomic_symbols[iz];
        // Case-insensitive comparison
        if (s.size() == sym.size()) {
            bool match = true;
            for (size_t i = 0; i < s.size(); ++i) {
                if (std::tolower(static_cast<unsigned char>(s[i])) !=
                    std::tolower(static_cast<unsigned char>(sym[i]))) {
                    match = false;
                    break;
                }
            }
            if (match) return iz;
        }
    }
    return 0;
}

} // namespace feff::common
