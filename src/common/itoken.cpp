// Token identifier for FEFF input file keywords.
// Converted from: src/COMMON/itoken.f

#include "itoken.hpp"
#include "string_utils.hpp"
#include <unordered_map>

namespace feff::common {

// Token maps keyed by the first 4 characters (uppercase).
// Built once on first call using function-local statics.

static const std::unordered_map<std::string, int>& feff_tokens() {
    static const std::unordered_map<std::string, int> map = {
        {"ATOM",  1},
        {"HOLE",  2},
        {"OVER",  3},
        {"CONT",  4},
        {"EXCH",  5},
        {"ION ",  6},
        {"TITL",  7},
        {"FOLP",  8},
        {"RPAT",  9},
        {"RMAX",  9},   // alias for RPAT
        {"DEBY", 10},
        {"RMUL", 11},
        {"SS  ", 12},
        {"PRIN", 13},
        {"POTE", 14},
        {"NLEG", 15},
        {"CRIT", 16},
        {"NOGE", 17},
        {"IORD", 18},
        {"PCRI", 19},
        {"SIG2", 20},
        {"XANE", 21},
        {"CORR", 22},
        {"AFOL", 23},
        {"EXAF", 24},
        {"POLA", 25},
        {"ELLI", 26},
        {"RGRI", 27},
        {"RPHA", 28},
        {"NSTA", 29},
        {"NOHO", 30},
        {"SIG3", 31},
        {"JUMP", 32},
        {"MBCO", 33},
        {"SPIN", 34},
        {"EDGE", 35},
        {"SCF ", 36},
        {"FMS ", 37},
        {"LDOS", 38},
        {"INTE", 39},
        {"CFAV", 40},
        {"S02 ", 41},
        {"XES ", 42},
        {"DANE", 43},
        {"FPRI", 44},
        {"RSIG", 45},
        {"XNCD", 46},
        {"XMCD", 46},   // alias for XNCD
        {"MULT", 47},
        {"UNFR", 48},
        {"TDLD", 49},
        {"PMBS", 50},
        {"PLAS", 51},
        {"SO2C", 52},
        {"SELF", 53},
        {"SFSE", 54},
        {"RCON", 55},   // RCONV -> first 4 chars = RCON
        {"ELNE", 56},
        {"EXEL", 57},
        {"MAGI", 58},
        {"ABSO", 59},
        {"EGRI", 60},
        {"END ", -1},
    };
    return map;
}

static const std::unordered_map<std::string, int>& spring_tokens() {
    static const std::unordered_map<std::string, int> map = {
        {"STRE",  1},
        {"ANGL",  2},
        {"VDOS",  3},
        {"PRDO",  4},   // PRDOS -> first 4 chars = PRDO
        {"END ", -1},
    };
    return map;
}

int itoken(const std::string& word, const std::string& filename) {
    // Extract first 4 characters, pad with spaces, uppercase
    std::string w = word.substr(0, 4);
    while (w.size() < 4) w += ' ';
    w = to_upper(w);

    // Select token map based on filename
    const std::unordered_map<std::string, int>* map = nullptr;

    if (filename.size() >= 8 && filename.substr(0, 8) == "feff.inp") {
        map = &feff_tokens();
    } else if (filename.size() >= 10 && filename.substr(0, 10) == "spring.inp") {
        map = &spring_tokens();
    } else {
        return 0;
    }

    auto it = map->find(w);
    return (it != map->end()) ? it->second : 0;
}

} // namespace feff::common
