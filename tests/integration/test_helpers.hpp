#pragma once
// Integration test utilities for comparing C++ pipeline output against
// Fortran baseline files.
//
// The work directory structure mirrors the Fortran test tree:
//   tests/integration/work/
//     Copper/
//       baseline/noSCF/    ← copied from Fortran baseline
//       baseline/withSCF/  ← copied from Fortran baseline
//       feff.inp           ← copied from baseline/noSCF/feff.inp
//       chi.dat            ← C++ pipeline output (when run)
//       xmu.dat            ← C++ pipeline output
//       feffNNNN.dat       ← C++ pipeline output
//     NiO/
//       ...

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#ifndef FORTRAN_TEST_DIR
#define FORTRAN_TEST_DIR "C:/Research_Software/feff85exafs_cpp/feff85exafs-fortran/tests"
#endif

#ifndef INTEGRATION_WORK_DIR
#define INTEGRATION_WORK_DIR "C:/Research_Software/feff85exafs_cpp/feff85exafs-cpp/tests/integration/work"
#endif

namespace fs = std::filesystem;

namespace feff::test {

// ---------------------------------------------------------------------------
// Path helpers
// ---------------------------------------------------------------------------

inline std::string get_fortran_test_dir() {
    return FORTRAN_TEST_DIR;
}

inline std::string get_work_dir() {
    return INTEGRATION_WORK_DIR;
}

/// Return the Fortran baseline directory for a material.
inline std::string get_fortran_baseline_dir(const std::string& material, bool scf = false) {
    return get_fortran_test_dir() + "/" + material + "/baseline/"
           + (scf ? "withSCF" : "noSCF");
}

/// Return the C++ work directory baseline for a material (mirrors Fortran structure).
inline std::string get_baseline_dir(const std::string& material, bool scf = false) {
    return get_work_dir() + "/" + material + "/baseline/"
           + (scf ? "withSCF" : "noSCF");
}

/// Return the C++ work directory root for a material (where pipeline runs).
inline std::string get_material_work_dir(const std::string& material) {
    return get_work_dir() + "/" + material;
}

inline bool file_exists(const std::string& path) {
    return fs::exists(path);
}

// ---------------------------------------------------------------------------
// Setup: mirror Fortran test structure into work directory
// ---------------------------------------------------------------------------

/// Set up the work directory for a material by copying the Fortran test
/// structure: baseline/noSCF/, baseline/withSCF/, and feff.inp.
inline void setup_material_work_dir(const std::string& material) {
    std::string fortran_material = get_fortran_test_dir() + "/" + material;
    std::string work_material = get_material_work_dir(material);

    // Create work directory
    fs::create_directories(work_material);

    // Copy baseline directories from Fortran
    for (const auto& scf_mode : {"noSCF", "withSCF"}) {
        std::string src_baseline = fortran_material + "/baseline/" + scf_mode;
        std::string dst_baseline = work_material + "/baseline/" + scf_mode;

        if (fs::is_directory(src_baseline)) {
            fs::create_directories(dst_baseline);
            // Copy all files from Fortran baseline
            for (const auto& entry : fs::directory_iterator(src_baseline)) {
                if (entry.is_regular_file()) {
                    fs::copy_file(entry.path(),
                                  fs::path(dst_baseline) / entry.path().filename(),
                                  fs::copy_options::overwrite_existing);
                }
            }
        }
    }

    // Copy feff.inp from baseline/noSCF/ to material root (for pipeline to use)
    std::string feff_src = work_material + "/baseline/noSCF/feff.inp";
    std::string feff_dst = work_material + "/feff.inp";
    if (fs::exists(feff_src)) {
        fs::copy_file(feff_src, feff_dst, fs::copy_options::overwrite_existing);
    }
}

/// Legacy helper — copies feff.inp into a given directory.
inline void copy_feff_inp(const std::string& material, const std::string& work_dir) {
    std::string src = get_fortran_baseline_dir(material) + "/feff.inp";
    std::string dst = work_dir + "/feff.inp";

    std::ifstream in(src, std::ios::binary);
    ASSERT_TRUE(in.good()) << "Cannot open source feff.inp: " << src;

    std::ofstream out(dst, std::ios::binary);
    ASSERT_TRUE(out.good()) << "Cannot create destination feff.inp: " << dst;

    out << in.rdbuf();
}

// ---------------------------------------------------------------------------
// Data file I/O
// ---------------------------------------------------------------------------

inline bool is_comment_line(const std::string& line) {
    auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return true;
    if (line[first] == '#') return true;
    std::istringstream iss(line.substr(first));
    double val;
    if (!(iss >> val)) return true;
    return false;
}

inline std::vector<std::vector<double>> read_data_file(const std::string& path,
                                                        int skip_header = 0) {
    std::vector<std::vector<double>> data;
    std::ifstream f(path);
    if (!f.good()) return data;

    std::string line;
    int skipped = 0;
    while (std::getline(f, line)) {
        if (is_comment_line(line)) continue;
        if (skipped < skip_header) { ++skipped; continue; }

        std::istringstream iss(line);
        std::vector<double> row;
        double val;
        while (iss >> val) row.push_back(val);
        if (!row.empty()) data.push_back(std::move(row));
    }
    return data;
}

// ---------------------------------------------------------------------------
// Comparison metrics
// ---------------------------------------------------------------------------

inline double r_factor(const std::vector<double>& ref,
                       const std::vector<double>& test) {
    if (ref.empty() || test.empty()) return 1.0e30;
    size_t n = std::min(ref.size(), test.size());
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = ref[i] - test[i];
        num += diff * diff;
        den += ref[i] * ref[i];
    }
    if (den < 1.0e-30) return 1.0e30;
    return num / den;
}

inline std::vector<double> extract_column(
        const std::vector<std::vector<double>>& data, int col) {
    std::vector<double> result;
    result.reserve(data.size());
    for (auto& row : data) {
        if (col < static_cast<int>(row.size()))
            result.push_back(row[col]);
    }
    return result;
}

inline void compare_data_files(const std::string& ref_file,
                               const std::string& test_file,
                               int col, double tol,
                               const std::string& label) {
    SCOPED_TRACE(label);
    ASSERT_TRUE(file_exists(ref_file)) << "Reference file missing: " << ref_file;
    ASSERT_TRUE(file_exists(test_file)) << "Test output file missing: " << test_file;

    auto ref_data = read_data_file(ref_file);
    auto test_data = read_data_file(test_file);

    ASSERT_FALSE(ref_data.empty()) << "No data rows in: " << ref_file;
    ASSERT_FALSE(test_data.empty()) << "No data rows in: " << test_file;

    auto ref_col = extract_column(ref_data, col);
    auto test_col = extract_column(test_data, col);

    double rf = r_factor(ref_col, test_col);
    EXPECT_LT(rf, tol)
        << label << ": R-factor " << rf << " exceeds " << tol
        << " (ref=" << ref_col.size() << " rows, test=" << test_col.size() << " rows)";
}

} // namespace feff::test
