Pipeline Module
===============

**Namespace:** ``feff::pipeline``

**Header:** ``src/pipeline/feff_pipeline.hpp``

The pipeline module orchestrates the complete FEFF calculation. It replaces the
Fortran main program ``feff.f`` and the individual ``ffmodN`` drivers.

PipelineConfig
--------------

.. code-block:: cpp

   struct PipelineConfig {
       bool run_pot = true;       // Stage 1+2: parse input + potentials
       bool run_xsph = true;      // Stage 3: phase shifts
       bool run_path = true;      // Stage 4: path finding
       bool run_genfmt = true;    // Stage 5: F-matrix
       bool run_ff2x = true;      // Stage 6: output spectra
       bool verbose = true;       // Enable progress logging
       std::string input_file = "feff.inp";
   };

Each boolean flag independently controls whether the corresponding pipeline
stage runs. Setting a flag to ``false`` skips only that stage, assuming
the required intermediate files already exist on disk. Note that RDINP (input
parsing) always runs regardless of flag settings.

Functions
---------

``int run_feff(const PipelineConfig& config)``
   Run the complete FEFF pipeline according to ``config``. Returns 0 on
   success, non-zero on error. This is the main entry point.

``int run_rdinp(const std::string& input_file)``
   Stage 1: Parse ``feff.inp`` and write JSON configuration files.

``int run_pot()``
   Stage 2: Calculate atomic potentials and write ``pot.pad``.

``int run_xsph()``
   Stage 3: Calculate phase shifts and write ``phase.pad``.

``int run_path()``
   Stage 4: Enumerate scattering paths and write ``paths.dat``, ``list.dat``.

``int run_genfmt()``
   Stage 5: Calculate F-matrices and write ``feff.pad``, ``feffNNNN.dat``.

``int run_ff2x()``
   Stage 6: Sum path contributions and write ``chi.dat``, ``xmu.dat``.

Example Usage
-------------

.. code-block:: cpp

   #include "src/pipeline/feff_pipeline.hpp"

   int main() {
       feff::pipeline::PipelineConfig config;
       config.input_file = "feff.inp";
       config.verbose = true;

       int result = feff::pipeline::run_feff(config);
       return result;
   }
