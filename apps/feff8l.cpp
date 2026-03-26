// feff8l — standalone FEFF8L EXAFS calculation executable.
//
// Usage: feff8l [feff.inp]
//
// Reads feff.inp (or specified input file), runs the complete FEFF
// calculation pipeline, and produces output files (chi.dat, xmu.dat,
// feffNNNN.dat, etc.)

#include "../src/pipeline/feff_pipeline.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    feff::pipeline::PipelineConfig config;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--skip-pot") {
            config.run_pot = false;
        } else if (arg == "--skip-to-path") {
            config.run_pot = false;
            config.run_xsph = false;
        } else if (arg == "--skip-to-genfmt") {
            config.run_pot = false;
            config.run_xsph = false;
            config.run_path = false;
        } else if (arg == "--skip-to-ff2x") {
            config.run_pot = false;
            config.run_xsph = false;
            config.run_path = false;
            config.run_genfmt = false;
        } else {
            config.input_file = arg;
        }
    }

    std::cout << "================================================" << std::endl;
    std::cout << " FEFF 8.5.0 (C++ port)" << std::endl;
    std::cout << " X-ray Absorption Spectroscopy Calculations" << std::endl;
    std::cout << "================================================" << std::endl;
    std::cout << std::endl;

    int result = feff::pipeline::run_feff(config);

    if (result != 0) {
        std::cerr << "FEFF exited with error code " << result << std::endl;
    }

    return result;
}
