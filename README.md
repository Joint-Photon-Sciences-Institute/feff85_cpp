# FEFF85 EXAFS - C++ Implementation

> **Status:** This C++ port is under active development and testing. It has been validated for **K-edge** and **L3-edge EXAFS** calculations, including **polarized** calculations. Testing against the original Fortran baseline is ongoing. **Use at your own risk** and always verify results against known references for your system.

A complete C++ port of the FEFF85 code for ab initio calculations of X-ray Absorption Fine Structure (EXAFS). This implementation reproduces the results of the original Fortran FEFF85 code and includes integration tests against the Fortran baseline for 8 materials.

## Prerequisites

**To build from source:**
- CMake 3.20 or later
- C++17 compiler (GCC/MinGW recommended on Windows, GCC/Clang on Linux/macOS)
- Ninja or Make build system

**To run the benchmark tests:**
- Python 3 with `numpy` and `matplotlib`

## Building from Source

```bash
# Clone the repository
git clone <repository-url>
cd feff85exafs-cpp-production

# Configure and build (Release mode)
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

# The executable will be at:
#   build/apps/feff8l.exe    (Windows)
#   build/apps/feff8l        (Linux/macOS)
```

The build statically links the C++ runtime, so the resulting executable is fully self-contained with no external DLL dependencies.

## Using Pre-built Binaries

If a pre-built binary is provided, simply place `feff8l` (or `feff8l.exe`) somewhere on your PATH. No additional libraries are required.

## Running FEFF

FEFF reads its input from a file named `feff.inp` in the current working directory:

```bash
cd /path/to/your/calculation
feff8l
```

The program will read `feff.inp` and produce output files including:
- `chi.dat` - EXAFS chi(k) data
- `xmu.dat` - X-ray absorption coefficient mu(E)
- `feffNNNN.dat` - Individual scattering path contributions
- `paths.dat` - Summary of scattering paths

## Benchmark Tests

The repository includes integration tests comparing the C++ output against the original Fortran FEFF85 baseline for 8 materials:

| Material | Description |
|---|---|
| Copper | FCC metal (Cu K-edge) |
| NiO | Nickel oxide (Ni K-edge) |
| UO2 | Uranium dioxide (U L3-edge) |
| LCO-para | La2CuO4, parallel polarization (Cu K-edge) |
| LCO-perp | La2CuO4, perpendicular polarization (Cu K-edge) |
| Zircon | ZrSiO4 (Zr K-edge) |
| bromoadamantane | C10H15Br (Br K-edge) |
| ferrocene | Fe(C5H5)2 (Fe K-edge) |

Each material is tested in both **SCF** (self-consistent field) and **noSCF** modes.

### Running the Benchmarks

```bash
# Build first (see above), then:
python tests/integration/benchmark.py
```

The script will:
1. Run `feff8l` in each material's `SCF/` and `noSCF/` folders
2. Compare the resulting `chi.dat` against the Fortran baseline
3. Compute R-factors for each comparison
4. Generate summary plots: `tests/integration/ALL_SCF_final.png` and `tests/integration/ALL_noSCF_final.png`
5. Print a pass/fail summary table (PASS if R-factor < 0.01)

## Project Structure

```
feff85exafs-cpp-production/
├── apps/                   # Standalone executable (feff8l)
├── cmake/                  # CMake helper modules
├── include/feff/           # Public API headers
├── src/                    # Source code (17 modules)
│   ├── atom/               # Atomic structure setup
│   ├── common/             # Shared constants and utilities
│   ├── debye/              # Debye-Waller factors
│   ├── exch/               # Exchange-correlation potentials
│   ├── ff2x/               # chi(k) and xmu(E) output
│   ├── fms/                # Full multiple scattering
│   ├── fovrg/              # Dirac equation solver
│   ├── genfmt/             # Scattering amplitude generation
│   ├── json_io/            # JSON I/O for intermediate data
│   ├── math/               # Math utilities (Bessel, matrices)
│   ├── par/                # Global parameters
│   ├── path/               # Path enumeration
│   ├── pipeline/           # Full calculation pipeline
│   ├── pot/                # Muffin-tin potentials
│   ├── rdinp/              # Input file parser (feff.inp)
│   └── xsph/               # Cross-section and phase shifts
├── tests/
│   ├── integration/        # Integration tests and benchmark script
│   │   ├── benchmark.py    # Benchmark runner
│   │   └── work/           # Test materials (8 materials)
│   │       └── <material>/
│   │           ├── SCF/        # SCF test input (feff.inp)
│   │           ├── noSCF/      # noSCF test input (feff.inp)
│   │           └── baseline/   # Fortran reference results
│   └── unit/               # Unit tests (Google Test)
└── third_party/            # External dependencies (Eigen, GoogleTest, nlohmann/json)
```

## License

See LICENSE file for details.
