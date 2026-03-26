#pragma once

// JSON I/O module for FEFF
// Converted from Fortran: RDINP/wrtjsn.f (writers) and JSON/*.f (readers)
// Uses nlohmann/json for JSON serialization
// JSON keys match Fortran exactly for cross-language compatibility

#include <array>
#include <complex>
#include <string>
#include <vector>

#include <feff/dimensions.hpp>
#include <feff/feff_input.hpp>

namespace feff::json_io {

// ============================================================================
// Writers (called from rdinp after parsing)
// Each writes a JSON file matching the Fortran output exactly
// ============================================================================

/// Write pot.json (MOD1 parameters: potentials, SCF, overlap)
void write_pot_json(const FeffInput& inp);

/// Write xsph.json (MOD2 parameters: phases, exchange)
void write_xsph_json(const FeffInput& inp);

/// Write path.json (MOD4 parameters: pathfinder)
void write_path_json(const FeffInput& inp);

/// Write genfmt.json (MOD5 parameters: F-matrix)
void write_genfmt_json(const FeffInput& inp);

/// Write ff2x.json (MOD6 parameters: chi output)
void write_ff2x_json(const FeffInput& inp);

/// Write atoms.json (atomic positions and potential types)
void write_atoms_json(const FeffInput& inp);

/// Write global.json (polarization, spin, configuration average)
void write_global_json(const FeffInput& inp, int nabs);

/// Write geom.json (geometry for pathfinder)
/// nat: number of atoms; rat[i] = {x,y,z}; iphat[i] = pot type; iatph[iph] = representative atom
void write_geom_json(const FeffInput& inp, int nat,
                     const double rat[][3], const int iphat[], const int iatph[]);

/// Write libpotph.json (combined input for libpotph library)
void write_libpotph_json(const FeffInput& inp);

/// Write xsect.json (cross-section data from XSPH)
void write_xsect_json(int ntit, const std::string titles[], double s02,
                       double erelax, double wp, double edge, double emu,
                       double gamach, int ne, int ne1, int ik0,
                       const double er[], const double ei[],
                       const double xsn[], const double col4[],
                       const double col5[],
                       const std::string& vfeff = vfeff_default,
                       const std::string& vf85e = vf85e_default);

/// Write feffNNNN.json (scattering path data from GENFMT)
/// fjson: output filename (e.g. "feff0001.json")
void write_feff_json(const std::string& fjson,
                     int ntit, const std::string titles[],
                     const double rat[][3], const int ipot[],
                     const double ri[], const double beta[], const double eta[],
                     int index, int iorder, int nleg, double deg,
                     double reff, double rnrmav, double edge, int ne,
                     const double col1[], const double col2[],
                     const double col3[], const double col4[],
                     const double col5[], const double col6[],
                     const double col7[],
                     const std::string& vfeff = vfeff_default,
                     const std::string& vf85e = vf85e_default);

// ============================================================================
// Readers (called by POT, XSPH, PATH, etc.)
// ============================================================================

/// Read atoms.json -> nat, rat(3,nat), iphat(nat)
void read_atoms_json(int& nat, double rat[][3], int iphat[]);

/// Read global.json -> all global parameters
void read_global_json(int& nabs, int& iphabs, double& rclabs,
                      int& ipol, int& ispin, int& le2,
                      double& elpty, double& angks,
                      double evec[3], double xivec[3], double spvec[3],
                      std::array<std::array<std::complex<double>, 3>, 3>& ptz);

/// Read geom.json -> geometry for pathfinder
void read_geom_json(int& nat, int& nph, int iatph[],
                    double rat[][3], int iphat[], int ibounc[]);

/// Read pot.json -> all MOD1 parameters
void read_pot_json(int& mpot, int& nph, int& ntitle, int& ihole,
                   int& ipr1, int& iafolp, int& ixc, int& ispec,
                   int& nmix, int& nohole, int& jumprm, int& inters,
                   int& nscmt, int& icoul, int& lfms1, int& iunf,
                   double& gamach, double& rgrd, double& ca1, double& ecv,
                   double& totvol, float& rfms1,
                   std::string title[], int iz[], int lmaxsc[],
                   double xnatph[], double xion[], double folp[],
                   int novr[], int iphovr[][novrx], int nnovr[][novrx],
                   double rovr[][novrx]);

/// Read xsect.json -> cross-section data
void read_xsect_json(int& ntit, std::string titles[],
                     double& s02, double& erelax, double& wp,
                     double& edge, double& emu, double& gamach,
                     int& ne, int& ne1, int& ik0,
                     double er[], double ei[], double xsn[],
                     double col4[], double col5[]);

/// Read only titles from xsect.json
void read_titles_json(int& ntit, std::string titles[]);

/// Read libpotph.json -> combined potentials + phases input
void read_libpotph_json(
    // TITLE
    int& ntitle, std::string title[],
    // ATOMS
    int& nat, double rat[][3], int iphat[],
    // POTENTIALS
    int& nph, int iz[], std::string potlbl[], int lmaxsc[], int lmaxph[],
    double xnatph[], double spinph[],
    // HOLE/EDGE
    int& ihole,
    // SCF
    float& rfms1, int& lfms1, int& nscmt, double& ca1, int& nmix,
    double& ecv, int& icoul,
    // POLARIZATION, ELLIPTICITY
    int& ipol, double evec[3], double& elpty, double xivec[3],
    // SPIN
    int& ispin, double spvec[3], double& angks,
    // computed
    std::array<std::array<std::complex<double>, 3>, 3>& ptz, double& gamach,
    // EXCHANGE
    int& ixc, double& vr0, double& vi0, int& ixc0,
    // AFOLP, FOLP, ION, RGRID, UNFREEZEF
    int& iafolp, double folp[], double xion[], double& rgrd, int& iunf,
    // INTERSTITIAL, JUMPRM, NOHOLE, PLASMON
    int& inters, double& totvol, int& jumprm, int& nohole, int& iplsmn);

// ============================================================================
// Utility
// ============================================================================

/// Throw a runtime_error if a required JSON key is missing
[[noreturn]] void bailout(const std::string& key, const std::string& file);

} // namespace feff::json_io
