// Read user-defined energy grid from grid.inp.
// Converted from src/XSPH/rdgrid.f

#include "rdgrid.hpp"
#include <feff/constants.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace feff::xsph {

// Helper: check if a string represents a number
static bool is_number(const std::string& s) {
    if (s.empty()) return false;
    size_t start = 0;
    if (s[0] == '+' || s[0] == '-') start = 1;
    if (start >= s.size()) return false;
    bool has_dot = false;
    bool has_e = false;
    for (size_t i = start; i < s.size(); i++) {
        if (std::isdigit(static_cast<unsigned char>(s[i]))) continue;
        if (s[i] == '.' && !has_dot) { has_dot = true; continue; }
        if ((s[i] == 'e' || s[i] == 'E' || s[i] == 'd' || s[i] == 'D') && !has_e) {
            has_e = true;
            if (i + 1 < s.size() && (s[i+1] == '+' || s[i+1] == '-')) i++;
            continue;
        }
        return false;
    }
    return true;
}

// Helper: skip comment lines
static bool is_comment(const std::string& line) {
    if (line.empty()) return true;
    char c = line[0];
    return (c == '#' || c == '!' || c == '*' || c == 'C' || c == 'c');
}

static void set_grid_min(double gridMin[], double gridMax[], double gridStep[],
                        int iGridType[], int nGrid) {
    int prev = nGrid - 1;
    int curr = nGrid;
    if ((iGridType[curr] != 2 && iGridType[prev] != 2) ||
        (iGridType[curr] == iGridType[prev])) {
        gridMin[curr] = gridMax[prev] + gridStep[curr];
    } else if (iGridType[curr] == 2) {
        gridMin[curr] = std::sqrt(2.0 * gridMax[prev] / hart) / bohr + gridStep[curr];
    } else {
        gridMin[curr] = (gridMax[prev] * bohr) * (gridMax[prev] * bohr) / 2.0 * hart + gridStep[curr];
    }
}

void rdgrid(FeffComplex em[], int& ne, int& nGrid,
            int iGridType[], double gridMin[], double gridMax[], double gridStep[],
            int nGridMax, int nex_dim) {

    std::ifstream fin("grid.inp");
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open grid.inp");
    }

    nGrid = 0;
    std::string line;

    while (nGrid < nGridMax) {
        // Skip comments
        while (std::getline(fin, line)) {
            if (!is_comment(line)) break;
        }
        if (fin.eof()) break;

        // Parse the line into words
        std::istringstream iss(line);
        std::string word1, word2, word3, word4;
        iss >> word1;
        if (word1.empty()) break;

        // Determine grid type
        std::string lower_w1 = word1;
        std::transform(lower_w1.begin(), lower_w1.end(), lower_w1.begin(), ::tolower);
        if (lower_w1 == "usergrid") {
            iGridType[nGrid] = 0;
        } else if (lower_w1 == "egrid") {
            iGridType[nGrid] = 1;
        } else if (lower_w1 == "kgrid") {
            iGridType[nGrid] = 2;
        } else if (lower_w1 == "expgrid") {
            iGridType[nGrid] = 3;
        } else {
            break;
        }

        if (iGridType[nGrid] != 0) {
            iss >> word2 >> word3 >> word4;
            if (word2 == "last") {
                if (nGrid > 0) {
                    set_grid_min(gridMin, gridMax, gridStep, iGridType, nGrid);
                } else {
                    gridMin[0] = 0.0;
                }
            } else {
                gridMin[nGrid] = std::stod(word2);
            }
            gridMax[nGrid] = std::stod(word3);
            gridStep[nGrid] = std::stod(word4);
        }

        if (iGridType[nGrid] == 0) {
            // User defined points
            while (std::getline(fin, line)) {
                if (is_comment(line)) continue;
                std::istringstream iss2(line);
                std::string w1, w2;
                iss2 >> w1;
                if (!is_number(w1)) {
                    // Back up for next grid type
                    gridMin[nGrid] = std::real(em[ne > 0 ? ne - 1 : 0]);
                    gridMax[nGrid] = std::real(em[ne > 0 ? ne - 1 : 0]);
                    break;
                }
                double realE = std::stod(w1);
                double imagE = 0.0;
                if (iss2 >> w2 && is_number(w2)) {
                    imagE = std::stod(w2);
                }
                em[ne] = FeffComplex(realE, imagE);
                ne++;
            }
        }
        nGrid++;
    }

    // Convert units
    for (int i = 0; i < nGrid; i++) {
        if (iGridType[i] == 2) {
            // k-grid: convert to inverse Bohr
            gridMin[i] *= bohr;
            gridMax[i] *= bohr;
            gridStep[i] *= bohr;
        } else if (iGridType[i] != 0) {
            // e-grid: convert to Hartrees
            gridMin[i] /= hart;
            gridMax[i] /= hart;
            gridStep[i] /= hart;
        }
    }
    // Convert user energies to Hartrees
    for (int i = 0; i < ne; i++) {
        em[i] = em[i] / hart;
    }
}

} // namespace feff::xsph
