// JSON I/O module for FEFF
// Converted from Fortran: RDINP/wrtjsn.f (writers) and JSON/*.f (readers)
// All JSON keys match the Fortran code exactly for cross-language compatibility.

#include <feff/json_io.hpp>
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/feff_input.hpp>

#include <nlohmann/json.hpp>

#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace feff::json_io {

// ============================================================================
// Utility
// ============================================================================

[[noreturn]] void bailout(const std::string& key, const std::string& file) {
    throw std::runtime_error(
        "Error reading " + file + ", failed to find >" + key + "<");
}

// Helper: get a value from JSON or call bailout
template <typename T>
static T jget(const json& j, const std::string& key, const std::string& file) {
    if (!j.contains(key)) bailout(key, file);
    return j.at(key).get<T>();
}

// Helper: write JSON to file with 2-space indent
static void write_json_file(const std::string& filename, const json& j) {
    std::ofstream ofs(filename);
    if (!ofs) {
        throw std::runtime_error("Cannot open " + filename + " for writing");
    }
    ofs << j.dump(2);
}

// Helper: read JSON from file
static json read_json_file(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("failed to read " + filename);
    }
    return json::parse(ifs);
}

// Helper: convert 0-based C++ array (size nphx+1) to JSON array
// Fortran writes iz(0:nphx) as a JSON array of length nphx+1
template <typename T, std::size_t N>
static std::vector<T> arr_to_vec(const std::array<T, N>& a) {
    return std::vector<T>(a.begin(), a.end());
}

template <typename T, std::size_t N>
static std::vector<T> arr_to_vec(const std::array<T, N>& a, int count) {
    return std::vector<T>(a.begin(), a.begin() + count);
}

// ============================================================================
// Writers
// ============================================================================

// --- pot.json (json_mod1) ---
void write_pot_json(const FeffInput& inp) {
    json j;
    j["mod"]    = 1;
    j["vfeff"]  = inp.vfeff;
    j["vf85e"]  = inp.vf85e;

    j["mpot"]   = inp.mpot;
    j["nph"]    = inp.nph;
    j["ntitle"] = inp.ntitle;
    j["ihole"]  = inp.ihole;
    j["ipr1"]   = inp.ipr1;
    j["iafolp"] = inp.iafolp;
    j["ixc"]    = inp.ixc;
    j["ispec"]  = inp.ispec;

    j["nmix"]   = inp.nmix;
    j["nohole"] = inp.nohole;
    j["jumprm"] = inp.jumprm;
    j["inters"] = inp.inters;
    j["nscmt"]  = inp.nscmt;
    j["icoul"]  = inp.icoul;
    j["lfms1"]  = inp.lfms1;
    j["iunf"]   = inp.iunf;

    j["gamach"] = inp.gamach;
    j["rgrd"]   = inp.rgrd;
    j["ca1"]    = inp.ca1;
    j["ecv"]    = inp.ecv;
    j["totvol"] = inp.totvol;
    j["rfms1"]  = static_cast<double>(inp.rfms1);

    // titles array (full nheadx)
    j["titles"] = arr_to_vec(inp.title);

    // 0:nphx arrays -> JSON arrays of length nphx+1
    j["iz"]     = arr_to_vec(inp.iz);
    j["lmaxsc"] = arr_to_vec(inp.lmaxsc);
    j["xnatph"] = arr_to_vec(inp.xnatph);
    j["xion"]   = arr_to_vec(inp.xion);
    j["folp"]   = arr_to_vec(inp.folp);
    j["novr"]   = arr_to_vec(inp.novr);

    // Overlap 2D arrays: disentangled per potential type
    // Fortran writes "iphovrN", "nnovrN", "rovrN" for each iph 0..nph
    for (int iph = 0; iph <= inp.nph; ++iph) {
        std::string suffix = std::to_string(iph);

        // iphovr(novrx, iph) -> "iphovrN"
        std::vector<int> iphovr_vec(novrx);
        for (int iovr = 0; iovr < novrx; ++iovr)
            iphovr_vec[iovr] = inp.iphovr[iph][iovr];
        j["iphovr" + suffix] = iphovr_vec;

        // nnovr(novrx, iph) -> "nnovrN"
        std::vector<int> nnovr_vec(novrx);
        for (int iovr = 0; iovr < novrx; ++iovr)
            nnovr_vec[iovr] = inp.nnovr[iph][iovr];
        j["nnovr" + suffix] = nnovr_vec;

        // rovr(novrx, iph) -> "rovrN"
        std::vector<double> rovr_vec(novrx);
        for (int iovr = 0; iovr < novrx; ++iovr)
            rovr_vec[iovr] = inp.rovr[iph][iovr];
        j["rovr" + suffix] = rovr_vec;
    }

    write_json_file("pot.json", j);
}

// --- xsph.json (json_mod2) ---
void write_xsph_json(const FeffInput& inp) {
    json j;
    j["mod"]      = 2;
    j["vfeff"]    = inp.vfeff;
    j["vf85e"]    = inp.vf85e;

    j["mphase"]   = inp.mphase;
    j["ipr2"]     = inp.ipr2;
    j["ixc"]      = inp.ixc;
    j["ixc0"]     = inp.ixc0;
    j["ispec"]    = inp.ispec;
    j["lreal"]    = inp.lreal;
    j["lfms2"]    = inp.lfms2;
    j["nph"]      = inp.nph;
    j["l2lp"]     = inp.l2lp;
    j["iPlsmn"]   = inp.iPlsmn;
    j["iGrid"]    = inp.iGrid;

    // Note: Fortran writes vr0/vi0 as keys "vro"/"vio"
    j["vro"]      = inp.vr0;
    j["vio"]      = inp.vi0;

    j["rgrd"]     = inp.rgrd;
    j["rfms2"]    = inp.rfms2;
    j["gamach"]   = inp.gamach;
    j["xkstep"]   = inp.xkstep;
    j["xkmax"]    = inp.xkmax;
    j["vixan"]    = inp.vixan;
    j["izstd"]    = inp.izstd;
    j["ifxc"]     = inp.ifxc;
    j["ipmbse"]   = inp.ipmbse;
    j["itdlda"]   = inp.itdlda;
    j["nonlocal"] = inp.nonlocal;
    j["ibasis"]   = inp.ibasis;

    j["lmaxph"]   = arr_to_vec(inp.lmaxph);

    // potlbl: array of strings
    std::vector<std::string> plbl(inp.potlbl.begin(), inp.potlbl.end());
    j["potlbl"]   = plbl;

    j["spinph"]   = arr_to_vec(inp.spinph);

    write_json_file("xsph.json", j);
}

// --- path.json (json_mod4) ---
void write_path_json(const FeffInput& inp) {
    json j;
    j["mod"]    = 4;
    j["vfeff"]  = inp.vfeff;
    j["vf85e"]  = inp.vf85e;

    j["mpath"]  = inp.mpath;
    j["ms"]     = inp.ms;
    j["nncrit"] = inp.nncrit;
    j["nlegxx"] = inp.nlegxx;
    j["ipr4"]   = inp.ipr4;
    j["critpw"] = static_cast<double>(inp.critpw);
    j["pcritk"] = static_cast<double>(inp.pcritk);
    j["pcrith"] = static_cast<double>(inp.pcrith);
    j["rmax"]   = static_cast<double>(inp.rmax);
    j["rfms2"]  = static_cast<double>(inp.rfms2);

    write_json_file("path.json", j);
}

// --- genfmt.json (json_mod5) ---
void write_genfmt_json(const FeffInput& inp) {
    json j;
    j["mod"]    = 5;
    j["vfeff"]  = inp.vfeff;
    j["vf85e"]  = inp.vf85e;

    j["mfeff"]  = inp.mfeff;
    j["ipr5"]   = inp.ipr5;
    j["iorder"] = inp.iorder;
    j["critcw"] = static_cast<double>(inp.critcw);
    j["wnstar"] = inp.wnstar;

    write_json_file("genfmt.json", j);
}

// --- ff2x.json (json_mod6) ---
void write_ff2x_json(const FeffInput& inp) {
    json j;
    j["mod"]    = 6;
    j["vfeff"]  = inp.vfeff;
    j["vf85e"]  = inp.vf85e;

    j["mchi"]   = inp.mchi;
    j["ispec"]  = inp.ispec;
    j["idwopt"] = inp.idwopt;
    j["ipr6"]   = inp.ipr6;
    j["mbconv"] = inp.mbconv;
    j["absolu"] = inp.absolu;
    j["vrcorr"] = inp.vrcorr;
    j["vicorr"] = inp.vicorr;
    j["s02"]    = inp.s02;
    j["critcw"] = inp.critcw;
    j["tk"]     = inp.tk;
    j["thetad"] = inp.thetad;
    j["alphat"] = inp.alphat;
    j["thetae"] = inp.thetae;
    j["sig2g"]  = inp.sig2g;

    write_json_file("ff2x.json", j);
}

// --- atoms.json (json_atoms) ---
void write_atoms_json(const FeffInput& inp) {
    json j;
    j["vfeff"] = inp.vfeff;
    j["vf85e"] = inp.vf85e;

    j["natt"] = inp.natt;

    std::vector<double> xx(inp.natt), yy(inp.natt), zz(inp.natt);
    std::vector<int> iph(inp.natt);
    for (int i = 0; i < inp.natt; ++i) {
        xx[i] = inp.ratx[i][0];
        yy[i] = inp.ratx[i][1];
        zz[i] = inp.ratx[i][2];
        iph[i] = inp.iphatx[i];
    }
    j["x"]      = xx;
    j["y"]      = yy;
    j["z"]      = zz;
    j["iphatx"] = iph;

    write_json_file("atoms.json", j);
}

// --- global.json (json_global) ---
void write_global_json(const FeffInput& inp, int nabs) {
    json j;
    j["vfeff"]  = inp.vfeff;
    j["vf85e"]  = inp.vf85e;

    j["nabs"]   = nabs;
    j["iphabs"] = inp.iphabs;
    j["rclabs"] = inp.rclabs;
    j["ipol"]   = inp.ipol;
    j["ispin"]  = inp.ispin;
    j["le2"]    = inp.le2;
    j["elpty"]  = inp.elpty;
    j["angks"]  = inp.angks;

    j["evec"]   = std::vector<double>{inp.evec[0], inp.evec[1], inp.evec[2]};
    j["xivec"]  = std::vector<double>{inp.xivec[0], inp.xivec[1], inp.xivec[2]};
    j["spvec"]  = std::vector<double>{inp.spvec[0], inp.spvec[1], inp.spvec[2]};

    // ptz(-1:1, -1:1) stored as ptz[row+1][col+1]
    // Fortran writes "ptz0" for i=-1 (i+1=0), "ptz1" for i=0, "ptz2" for i=1
    // Each ptzN is a 6-element array: [real(-1,i), imag(-1,i), real(0,i), imag(0,i), real(1,i), imag(1,i)]
    for (int i = -1; i <= 1; ++i) {
        std::string key = "ptz" + std::to_string(i + 1);
        std::vector<double> pv(6);
        pv[0] = inp.ptz[0][i + 1].real();  // ptz(-1, i)
        pv[1] = inp.ptz[0][i + 1].imag();
        pv[2] = inp.ptz[1][i + 1].real();  // ptz( 0, i)
        pv[3] = inp.ptz[1][i + 1].imag();
        pv[4] = inp.ptz[2][i + 1].real();  // ptz( 1, i)
        pv[5] = inp.ptz[2][i + 1].imag();
        j[key] = pv;
    }

    write_json_file("global.json", j);
}

// --- geom.json (json_geom) ---
void write_geom_json(const FeffInput& inp, int nat,
                     const double rat[][3], const int iphat[], const int iatph[]) {
    json j;
    j["vfeff"] = inp.vfeff;
    j["vf85e"] = inp.vf85e;

    j["natt"] = inp.natt;
    j["nph"]  = inp.nph;

    // iatph(0:nphx) -> JSON array of length nphx+1
    std::vector<int> iatph_vec(nphx + 1);
    for (int iph = 0; iph <= nphx; ++iph)
        iatph_vec[iph] = iatph[iph];
    j["iatph"] = iatph_vec;

    // Atom coordinates and types (1-based in Fortran, 0-based here)
    std::vector<double> xx(inp.natt), yy(inp.natt), zz(inp.natt);
    std::vector<int> iph_vec(inp.natt), ibo_vec(inp.natt);
    for (int i = 0; i < inp.natt; ++i) {
        xx[i] = rat[i][0];
        yy[i] = rat[i][1];
        zz[i] = rat[i][2];
        iph_vec[i] = iphat[i];
        ibo_vec[i] = 1;  // Fortran always sets ib=1
    }
    j["x"]   = xx;
    j["y"]   = yy;
    j["z"]   = zz;
    j["iph"] = iph_vec;
    j["ibo"] = ibo_vec;

    write_json_file("geom.json", j);
}

// --- libpotph.json (json_libpotph) ---
void write_libpotph_json(const FeffInput& inp) {
    json j;

    // Version
    j["vfeff"] = inp.vfeff;
    j["vf85e"] = inp.vf85e;

    // TITLE
    j["ntitle"] = inp.ntitle;
    j["titles"] = arr_to_vec(inp.title);

    // ATOMS
    j["natt"] = inp.natt;
    std::vector<double> xx(inp.natt), yy(inp.natt), zz(inp.natt);
    std::vector<int> iph(inp.natt);
    for (int i = 0; i < inp.natt; ++i) {
        xx[i] = inp.ratx[i][0];
        yy[i] = inp.ratx[i][1];
        zz[i] = inp.ratx[i][2];
        iph[i] = inp.iphatx[i];
    }
    j["x"]      = xx;
    j["y"]      = yy;
    j["z"]      = zz;
    j["iphatx"] = iph;

    // POTENTIALS
    j["nph"]    = inp.nph;
    j["iz"]     = arr_to_vec(inp.iz);
    std::vector<std::string> plbl(inp.potlbl.begin(), inp.potlbl.end());
    j["potlbl"] = plbl;
    j["lmaxsc"] = arr_to_vec(inp.lmaxsc);
    j["lmaxph"] = arr_to_vec(inp.lmaxph);
    j["xnatph"] = arr_to_vec(inp.xnatph);
    j["spinph"] = arr_to_vec(inp.spinph);

    // HOLE/EDGE
    j["ihole"] = inp.ihole;

    // SCF
    j["rfms1"] = static_cast<double>(inp.rfms1);
    j["lfms1"] = inp.lfms1;
    j["nscmt"] = inp.nscmt;
    j["ca1"]   = inp.ca1;
    j["nmix"]  = inp.nmix;
    j["ecv"]   = inp.ecv;
    j["icoul"] = inp.icoul;

    // POLARIZATION, ELLIPTICITY
    j["ipol"]  = inp.ipol;
    j["evec"]  = std::vector<double>{inp.evec[0], inp.evec[1], inp.evec[2]};
    j["elpty"] = inp.elpty;
    j["xivec"] = std::vector<double>{inp.xivec[0], inp.xivec[1], inp.xivec[2]};

    // SPIN
    j["ispin"] = inp.ispin;
    j["spvec"] = std::vector<double>{inp.spvec[0], inp.spvec[1], inp.spvec[2]};
    j["angks"] = inp.angks;

    // computed: ptz and gamach
    for (int i = -1; i <= 1; ++i) {
        std::string key = "ptz" + std::to_string(i + 1);
        std::vector<double> pv(6);
        pv[0] = inp.ptz[0][i + 1].real();
        pv[1] = inp.ptz[0][i + 1].imag();
        pv[2] = inp.ptz[1][i + 1].real();
        pv[3] = inp.ptz[1][i + 1].imag();
        pv[4] = inp.ptz[2][i + 1].real();
        pv[5] = inp.ptz[2][i + 1].imag();
        j[key] = pv;
    }
    j["gamach"] = inp.gamach;

    // EXCHANGE
    j["ixc"]  = inp.ixc;
    j["vro"]  = inp.vr0;   // Fortran key is "vro" not "vr0"
    j["vio"]  = inp.vi0;   // Fortran key is "vio" not "vi0"
    j["ixc0"] = inp.ixc0;

    // AFOLP, FOLP, ION, RGRID, UNFREEZEF
    j["iafolp"] = inp.iafolp;
    j["folp"]   = arr_to_vec(inp.folp);
    j["xion"]   = arr_to_vec(inp.xion);
    j["rgrd"]   = inp.rgrd;
    j["iunf"]   = inp.iunf;

    // INTERSTITIAL, JUMPRM, NOHOLE
    j["inters"] = inp.inters;
    j["totvol"] = inp.totvol;
    j["jumprm"] = inp.jumprm;
    j["nohole"] = inp.nohole;
    j["iplsmn"] = inp.iPlsmn;

    write_json_file("libpotph.json", j);
}

// --- xsect.json (json_xsect) ---
void write_xsect_json(int ntit, const std::string titles[], double s02,
                       double erelax, double wp, double edge, double emu,
                       double gamach, int ne, int ne1, int ik0,
                       const double er[], const double ei[],
                       const double xsn[], const double col4[],
                       const double col5[],
                       const std::string& vfeff_str,
                       const std::string& vf85e_str) {
    json j;
    j["vfeff"]  = vfeff_str;
    j["vf85e"]  = vf85e_str;

    j["ntitle"] = ntit;
    j["title"]  = std::vector<std::string>(titles, titles + ntit);

    j["s02"]    = s02;
    j["erelax"] = erelax;
    j["wp"]     = wp;
    j["edge"]   = edge;
    j["emu"]    = emu;
    j["gamach"] = gamach;
    j["ne"]     = ne;
    j["ne1"]    = ne1;
    j["ik0"]    = ik0;

    j["ereal"]  = std::vector<double>(er,   er   + ne);
    j["eimag"]  = std::vector<double>(ei,   ei   + ne);
    j["xsnorm"] = std::vector<double>(xsn,  xsn  + ne);
    j["dum1"]   = std::vector<double>(col4, col4 + ne);
    j["dum2"]   = std::vector<double>(col5, col5 + ne);

    write_json_file("xsect.json", j);
}

// --- feffNNNN.json (json_nnnn) ---
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
                     const std::string& vfeff_str,
                     const std::string& vf85e_str) {
    constexpr double eps = 1.0e-8;

    json j;
    j["vfeff"]   = vfeff_str;
    j["vf85e"]   = vf85e_str;
    j["titles"]  = std::vector<std::string>(titles, titles + ntit);

    j["index"]   = index;
    j["iorder"]  = iorder;
    j["nleg"]    = nleg;
    j["degen"]   = deg;
    j["reff"]    = reff * bohr;     // convert from Bohr to Angstrom
    j["rnorman"] = rnrmav;
    j["edge"]    = edge * hart;     // convert from Hartree to eV

    // Leg lengths, beta and eta angles (convert to Angstrom / degrees)
    std::vector<double> ria(nleg), betad(nleg), etad(nleg);
    for (int i = 0; i < nleg; ++i) {
        ria[i]   = ri[i] * bohr;
        betad[i] = beta[i] * 180.0 / pi;
        etad[i]  = eta[i] * 180.0 / pi;
        if (std::abs(etad[i] - 360.0) < eps) {
            etad[i] = 0.0;
        }
    }
    j["ri"]   = ria;
    j["beta"] = betad;
    j["eta"]  = etad;

    // Atom positions: "atomN" where N=1..nleg
    // Each is a 4-element array [x, y, z, ipot] with coords in Angstrom
    for (int iat = 0; iat < nleg; ++iat) {
        std::string key = "atom" + std::to_string(iat + 1);
        std::vector<double> atom(4);
        atom[0] = rat[iat][0] * bohr;
        atom[1] = rat[iat][1] * bohr;
        atom[2] = rat[iat][2] * bohr;
        atom[3] = static_cast<double>(ipot[iat]);
        j[key] = atom;
    }

    // Data columns
    j["k"]        = std::vector<double>(col1, col1 + ne);
    j["real_phc"] = std::vector<double>(col2, col2 + ne);
    j["mag_feff"] = std::vector<double>(col3, col3 + ne);
    j["pha_feff"] = std::vector<double>(col4, col4 + ne);
    j["red_fact"] = std::vector<double>(col5, col5 + ne);
    j["lam"]      = std::vector<double>(col6, col6 + ne);
    j["rep"]      = std::vector<double>(col7, col7 + ne);

    write_json_file(fjson, j);
}


// ============================================================================
// Readers
// ============================================================================

// --- read_atoms_json (json_read_atoms) ---
void read_atoms_json(int& nat, double rat[][3], int iphat[]) {
    const std::string file = "atoms.json";
    json j = read_json_file(file);

    nat = jget<int>(j, "natt", file);

    auto xx    = jget<std::vector<double>>(j, "x", file);
    auto yy    = jget<std::vector<double>>(j, "y", file);
    auto zz    = jget<std::vector<double>>(j, "z", file);
    auto iphvec = jget<std::vector<int>>(j, "iphatx", file);

    for (int i = 0; i < nat; ++i) {
        iphat[i]  = iphvec[i];
        rat[i][0] = xx[i];
        rat[i][1] = yy[i];
        rat[i][2] = zz[i];
    }
}

// --- read_global_json (json_read_global) ---
void read_global_json(int& nabs, int& iphabs, double& rclabs,
                      int& ipol, int& ispin, int& le2,
                      double& elpty, double& angks,
                      double evec[3], double xivec[3], double spvec[3],
                      std::array<std::array<std::complex<double>, 3>, 3>& ptz) {
    const std::string file = "global.json";
    json j = read_json_file(file);

    nabs   = jget<int>(j, "nabs", file);
    iphabs = jget<int>(j, "iphabs", file);
    rclabs = jget<double>(j, "rclabs", file);
    ipol   = jget<int>(j, "ipol", file);
    ispin  = jget<int>(j, "ispin", file);
    le2    = jget<int>(j, "le2", file);
    elpty  = jget<double>(j, "elpty", file);
    angks  = jget<double>(j, "angks", file);

    auto ev = jget<std::vector<double>>(j, "evec", file);
    for (int i = 0; i < 3; ++i) evec[i] = ev[i];

    auto xi = jget<std::vector<double>>(j, "xivec", file);
    for (int i = 0; i < 3; ++i) xivec[i] = xi[i];

    auto sp = jget<std::vector<double>>(j, "spvec", file);
    for (int i = 0; i < 3; ++i) spvec[i] = sp[i];

    // ptz: "ptz0" -> column i=-1 (index 0), "ptz1" -> column i=0 (index 1), etc.
    // Each ptzN has 6 doubles: [re(-1,i), im(-1,i), re(0,i), im(0,i), re(1,i), im(1,i)]
    for (int i = -1; i <= 1; ++i) {
        std::string key = "ptz" + std::to_string(i + 1);
        auto pv = jget<std::vector<double>>(j, key, file);
        ptz[0][i + 1] = std::complex<double>(pv[0], pv[1]);  // ptz(-1, i)
        ptz[1][i + 1] = std::complex<double>(pv[2], pv[3]);  // ptz( 0, i)
        ptz[2][i + 1] = std::complex<double>(pv[4], pv[5]);  // ptz( 1, i)
    }
}

// --- read_geom_json (json_read_geom) ---
void read_geom_json(int& nat, int& nph, int iatph[],
                    double rat[][3], int iphat[], int ibounc[]) {
    const std::string file = "geom.json";

    // Initialize iatph
    nph = 0;
    for (int iph = 0; iph <= nphx; ++iph)
        iatph[iph] = 0;

    json j = read_json_file(file);

    nat = jget<int>(j, "natt", file);

    auto xx  = jget<std::vector<double>>(j, "x", file);
    auto yy  = jget<std::vector<double>>(j, "y", file);
    auto zz  = jget<std::vector<double>>(j, "z", file);
    auto iph = jget<std::vector<int>>(j, "iph", file);
    auto ibo = jget<std::vector<int>>(j, "ibo", file);

    for (int i = 0; i < nat; ++i) {
        iphat[i]  = iph[i];
        rat[i][0] = xx[i];
        rat[i][1] = yy[i];
        rat[i][2] = zz[i];
        ibounc[i] = ibo[i];

        // Reconstruct nph and iatph like the Fortran code
        if (iphat[i] > nph) nph = iphat[i];
        if (iatph[iphat[i]] == 0) iatph[iphat[i]] = i + 1;  // 1-based atom index
    }
}

// --- read_pot_json (json_read_pot) ---
void read_pot_json(int& mpot, int& nph, int& ntitle, int& ihole,
                   int& ipr1, int& iafolp, int& ixc, int& ispec,
                   int& nmix, int& nohole, int& jumprm, int& inters,
                   int& nscmt, int& icoul, int& lfms1, int& iunf,
                   double& gamach, double& rgrd, double& ca1, double& ecv,
                   double& totvol, float& rfms1,
                   std::string title[], int iz[], int lmaxsc[],
                   double xnatph[], double xion[], double folp[],
                   int novr[], int iphovr[][novrx], int nnovr[][novrx],
                   double rovr[][novrx]) {
    const std::string file = "pot.json";
    json j = read_json_file(file);

    mpot   = jget<int>(j, "mpot", file);
    nph    = jget<int>(j, "nph", file);
    ntitle = jget<int>(j, "ntitle", file);
    ihole  = jget<int>(j, "ihole", file);
    ipr1   = jget<int>(j, "ipr1", file);
    iafolp = jget<int>(j, "iafolp", file);
    ixc    = jget<int>(j, "ixc", file);
    ispec  = jget<int>(j, "ispec", file);

    nmix   = jget<int>(j, "nmix", file);
    nohole = jget<int>(j, "nohole", file);
    jumprm = jget<int>(j, "jumprm", file);
    inters = jget<int>(j, "inters", file);
    nscmt  = jget<int>(j, "nscmt", file);
    icoul  = jget<int>(j, "icoul", file);
    lfms1  = jget<int>(j, "lfms1", file);
    iunf   = jget<int>(j, "iunf", file);

    gamach = jget<double>(j, "gamach", file);
    rgrd   = jget<double>(j, "rgrd", file);
    ca1    = jget<double>(j, "ca1", file);
    ecv    = jget<double>(j, "ecv", file);
    totvol = jget<double>(j, "totvol", file);
    rfms1  = static_cast<float>(jget<double>(j, "rfms1", file));

    // titles: array of nheadx strings
    auto strings = jget<std::vector<std::string>>(j, "titles", file);
    for (int i = 0; i < nheadx && i < static_cast<int>(strings.size()); ++i)
        title[i] = strings[i];

    // 0:nphx arrays (JSON stores as 0-indexed vectors of length nphx+1)
    auto iz_vec     = jget<std::vector<int>>(j, "iz", file);
    auto lmaxsc_vec = jget<std::vector<int>>(j, "lmaxsc", file);
    auto xnatph_vec = jget<std::vector<double>>(j, "xnatph", file);
    auto xion_vec   = jget<std::vector<double>>(j, "xion", file);
    auto folp_vec   = jget<std::vector<double>>(j, "folp", file);
    auto novr_vec   = jget<std::vector<int>>(j, "novr", file);

    for (int iph = 0; iph <= nphx; ++iph) {
        iz[iph]     = iz_vec[iph];
        lmaxsc[iph] = lmaxsc_vec[iph];
        xnatph[iph] = xnatph_vec[iph];
        xion[iph]   = xion_vec[iph];
        folp[iph]   = folp_vec[iph];
        novr[iph]   = novr_vec[iph];
    }

    // Overlap 2D arrays
    for (int iph = 0; iph <= nph; ++iph) {
        std::string suffix = std::to_string(iph);

        auto iphovr_vec = jget<std::vector<int>>(j, "iphovr" + suffix, file);
        auto nnovr_vec2 = jget<std::vector<int>>(j, "nnovr" + suffix, file);
        auto rovr_vec   = jget<std::vector<double>>(j, "rovr" + suffix, file);

        for (int iovr = 0; iovr < novr[iph]; ++iovr) {
            iphovr[iph][iovr] = iphovr_vec[iovr];
            nnovr[iph][iovr]  = nnovr_vec2[iovr];
            rovr[iph][iovr]   = rovr_vec[iovr];
        }
    }
}

// --- read_xsect_json (read_xsect) ---
void read_xsect_json(int& ntit, std::string titles[],
                     double& s02, double& erelax, double& wp,
                     double& edge, double& emu, double& gamach,
                     int& ne, int& ne1, int& ik0,
                     double er[], double ei[], double xsn[],
                     double col4[], double col5[]) {
    const std::string file = "xsect.json";
    json j = read_json_file(file);

    ntit = jget<int>(j, "ntitle", file);
    auto strings = jget<std::vector<std::string>>(j, "title", file);
    for (int i = 0; i < ntit; ++i)
        titles[i] = strings[i];

    s02    = jget<double>(j, "s02", file);
    erelax = jget<double>(j, "erelax", file);
    wp     = jget<double>(j, "wp", file);
    edge   = jget<double>(j, "edge", file);
    emu    = jget<double>(j, "emu", file);
    gamach = jget<double>(j, "gamach", file);
    ne     = jget<int>(j, "ne", file);
    ne1    = jget<int>(j, "ne1", file);
    ik0    = jget<int>(j, "ik0", file);

    auto er_vec  = jget<std::vector<double>>(j, "ereal", file);
    auto ei_vec  = jget<std::vector<double>>(j, "eimag", file);
    auto xsn_vec = jget<std::vector<double>>(j, "xsnorm", file);
    auto c4_vec  = jget<std::vector<double>>(j, "dum1", file);
    auto c5_vec  = jget<std::vector<double>>(j, "dum2", file);

    for (int i = 0; i < ne; ++i) {
        er[i]   = er_vec[i];
        ei[i]   = ei_vec[i];
        xsn[i]  = xsn_vec[i];
        col4[i] = c4_vec[i];
        col5[i] = c5_vec[i];
    }
}

// --- read_titles_json (read_titles) ---
void read_titles_json(int& ntit, std::string titles[]) {
    const std::string file = "xsect.json";
    json j = read_json_file(file);

    ntit = jget<int>(j, "ntitle", file);
    auto strings = jget<std::vector<std::string>>(j, "title", file);
    for (int i = 0; i < ntit; ++i)
        titles[i] = strings[i];
}

// --- read_libpotph_json (json_read_libpotph) ---
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
    int& inters, double& totvol, int& jumprm, int& nohole, int& iplsmn) {

    const std::string file = "libpotph.json";
    json j = read_json_file(file);

    // TITLE
    ntitle = jget<int>(j, "ntitle", file);
    auto strings = jget<std::vector<std::string>>(j, "titles", file);
    for (int i = 0; i < nheadx && i < static_cast<int>(strings.size()); ++i)
        title[i] = strings[i];

    // ATOMS
    nat = jget<int>(j, "natt", file);
    auto xx     = jget<std::vector<double>>(j, "x", file);
    auto yy     = jget<std::vector<double>>(j, "y", file);
    auto zz     = jget<std::vector<double>>(j, "z", file);
    auto ip_vec = jget<std::vector<int>>(j, "iphatx", file);
    for (int i = 0; i < nat; ++i) {
        iphat[i]  = ip_vec[i];
        rat[i][0] = xx[i];
        rat[i][1] = yy[i];
        rat[i][2] = zz[i];
    }

    // POTENTIALS
    nph = jget<int>(j, "nph", file);

    auto iz_vec     = jget<std::vector<int>>(j, "iz", file);
    for (int iph = 0; iph <= nphx; ++iph)
        iz[iph] = iz_vec[iph];

    auto plbl_vec = jget<std::vector<std::string>>(j, "potlbl", file);
    for (int iph = 0; iph <= nphx && iph < static_cast<int>(plbl_vec.size()); ++iph)
        potlbl[iph] = plbl_vec[iph].substr(0, 6);  // Fortran truncates to 6 chars

    auto lmaxsc_vec = jget<std::vector<int>>(j, "lmaxsc", file);
    for (int iph = 0; iph <= nphx; ++iph)
        lmaxsc[iph] = lmaxsc_vec[iph];

    auto lmaxph_vec = jget<std::vector<int>>(j, "lmaxph", file);
    for (int iph = 0; iph <= nphx; ++iph)
        lmaxph[iph] = lmaxph_vec[iph];

    auto xnatph_vec = jget<std::vector<double>>(j, "xnatph", file);
    for (int iph = 0; iph <= nphx; ++iph)
        xnatph[iph] = xnatph_vec[iph];

    auto spinph_vec = jget<std::vector<double>>(j, "spinph", file);
    for (int iph = 0; iph <= nphx; ++iph)
        spinph[iph] = spinph_vec[iph];

    // HOLE/EDGE
    ihole = jget<int>(j, "ihole", file);

    // SCF
    rfms1 = static_cast<float>(jget<double>(j, "rfms1", file));
    lfms1 = jget<int>(j, "lfms1", file);
    nscmt = jget<int>(j, "nscmt", file);
    ca1   = jget<double>(j, "ca1", file);
    nmix  = jget<int>(j, "nmix", file);
    ecv   = jget<double>(j, "ecv", file);
    icoul = jget<int>(j, "icoul", file);

    // POLARIZATION, ELLIPTICITY
    ipol = jget<int>(j, "ipol", file);
    auto ev = jget<std::vector<double>>(j, "evec", file);
    for (int i = 0; i < 3; ++i) evec[i] = ev[i];
    elpty = jget<double>(j, "elpty", file);
    auto xi = jget<std::vector<double>>(j, "xivec", file);
    for (int i = 0; i < 3; ++i) xivec[i] = xi[i];

    // SPIN
    ispin = jget<int>(j, "ispin", file);
    auto sp = jget<std::vector<double>>(j, "spvec", file);
    for (int i = 0; i < 3; ++i) spvec[i] = sp[i];
    angks = jget<double>(j, "angks", file);

    // computed: ptz and gamach
    for (int i = -1; i <= 1; ++i) {
        std::string key = "ptz" + std::to_string(i + 1);
        auto pv = jget<std::vector<double>>(j, key, file);
        ptz[0][i + 1] = std::complex<double>(pv[0], pv[1]);  // ptz(-1, i)
        ptz[1][i + 1] = std::complex<double>(pv[2], pv[3]);  // ptz( 0, i)
        ptz[2][i + 1] = std::complex<double>(pv[4], pv[5]);  // ptz( 1, i)
    }
    gamach = jget<double>(j, "gamach", file);

    // EXCHANGE
    ixc  = jget<int>(j, "ixc", file);
    vr0  = jget<double>(j, "vro", file);   // JSON key is "vro"
    vi0  = jget<double>(j, "vio", file);   // JSON key is "vio"
    ixc0 = jget<int>(j, "ixc0", file);

    // AFOLP, FOLP, ION, RGRID, UNFREEZEF
    iafolp = jget<int>(j, "iafolp", file);
    auto folp_vec = jget<std::vector<double>>(j, "folp", file);
    for (int iph = 0; iph <= nphx; ++iph)
        folp[iph] = folp_vec[iph];
    auto xion_vec = jget<std::vector<double>>(j, "xion", file);
    for (int iph = 0; iph <= nphx; ++iph)
        xion[iph] = xion_vec[iph];
    rgrd = jget<double>(j, "rgrd", file);
    iunf = jget<int>(j, "iunf", file);

    // INTERSTITIAL, JUMPRM, NOHOLE, PLASMON
    inters = jget<int>(j, "inters", file);
    totvol = jget<double>(j, "totvol", file);
    nohole = jget<int>(j, "nohole", file);
    jumprm = jget<int>(j, "jumprm", file);
    iplsmn = jget<int>(j, "iplsmn", file);

    // Unit conversions (matching Fortran json_read_libpotph lines 342-350)
    rfms1  = rfms1 / static_cast<float>(bohr);
    gamach = gamach / hart;
    ecv    = ecv / hart;
    totvol = totvol / (bohr * bohr * bohr);
    for (int iat = 0; iat < nat; ++iat) {
        rat[iat][0] /= bohr;
        rat[iat][1] /= bohr;
        rat[iat][2] /= bohr;
    }
}

} // namespace feff::json_io
