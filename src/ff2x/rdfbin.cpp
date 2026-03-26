// Read feff.pad (Packed ASCII Data format).
// Converted from: src/FF2X/rdfbin.f
//
// PAD format: printable ASCII with special line markers:
//   #_    header lines (version, npot, ne)
//   #"    title / plain text
//   #&    misc info (ihole, iorder, ilinit, rnrmav, xmu, edge)
//   #@    potential labels and iz
//   ##    path info (index, nleg, deg, reff, crit, ipots)
//   !     PAD-encoded real array
//   $     PAD-encoded complex array

#include "rdfbin.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "../common/logging.hpp"
#include "../common/pad_io.hpp"
#include "../common/string_utils.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace feff::ff2x {

FeffPadData read_feff_pad(const std::string& filename) {
    FeffPadData d;

    std::string fname = filename;
    if (fname.empty()) fname = "feff.pad";

    std::ifstream fin(fname);
    if (!fin.is_open()) {
        // Not an error for optional files -- caller decides
        return d;
    }

    std::string str;

    // First line: must start with "#_feff.pad"
    if (!std::getline(fin, str)) {
        throw std::runtime_error("rdfbin error: cannot read first line of " + fname);
    }
    // Trim leading whitespace
    auto pos = str.find_first_not_of(" \t");
    if (pos != std::string::npos) str = str.substr(pos);

    if (str.substr(0, 10) != "#_feff.pad") {
        throw std::runtime_error("rdfbin error: wrong format in " + fname);
    }

    // Check version -- only v03 supported
    int ivers = 0;
    if (str.size() >= 14) {
        std::string tag = str.substr(0, 14);
        if (tag == "#_feff.pad fil") ivers = 1;
        else if (tag == "#_feff.pad v02") ivers = 2;
        else if (tag == "#_feff.pad v03") ivers = 3;
    }
    if (ivers != 3) {
        throw std::runtime_error("rdfbin error: unsupported version of feff.pad");
    }

    // Second line: #_ npot ne npadx
    if (!std::getline(fin, str)) {
        throw std::runtime_error("rdfbin error: premature end of " + fname);
    }
    pos = str.find_first_not_of(" \t");
    if (pos != std::string::npos) str = str.substr(pos);
    if (str.substr(0, 2) != "#_") {
        throw std::runtime_error("rdfbin error: expected #_ line in " + fname);
    }

    int npot_read = 0, ne_read = 0, npadx = 0;
    {
        std::istringstream iss(str.substr(2));
        if (!(iss >> npot_read >> ne_read >> npadx)) {
            throw std::runtime_error("rdfbin error: bad data on line 2 of " + fname);
        }
    }
    d.npot = npot_read;
    d.ne = ne_read;

    // Third line: #& ihole iorder ilinit rnrmav xmu edge
    if (!std::getline(fin, str)) {
        throw std::runtime_error("rdfbin error: premature end of " + fname);
    }
    pos = str.find_first_not_of(" \t");
    if (pos != std::string::npos) str = str.substr(pos);
    if (str.substr(0, 2) != "#&") {
        throw std::runtime_error("rdfbin error: expected #& line in " + fname);
    }
    {
        std::istringstream iss(str.substr(2));
        if (!(iss >> d.ihole >> d.iorder >> d.ilinit >> d.rnrmav >> d.xmu >> d.edge)) {
            throw std::runtime_error("rdfbin error: bad data on #& line in " + fname);
        }
    }

    // Fourth line: #@ potlbl(0:npot) iz(0:npot)
    if (!std::getline(fin, str)) {
        throw std::runtime_error("rdfbin error: premature end of " + fname);
    }
    pos = str.find_first_not_of(" \t");
    if (pos != std::string::npos) str = str.substr(pos);
    if (str.substr(0, 2) != "#@") {
        throw std::runtime_error("rdfbin error: expected #@ line in " + fname);
    }
    {
        std::istringstream iss(str.substr(2));
        // Read potlbl(0:npot) then iz(0:npot)
        for (int i = 0; i <= d.npot; ++i) {
            if (!(iss >> d.potlbl[i])) {
                throw std::runtime_error("rdfbin error: missing pot label in " + fname);
            }
        }
        for (int i = 0; i <= d.npot; ++i) {
            if (!(iss >> d.iz[i])) {
                throw std::runtime_error("rdfbin error: missing iz in " + fname);
            }
        }
    }

    // Read PAD-encoded arrays common to all paths
    // phc(ne) -- complex
    {
        std::complex<float> tmp_c[nex];
        // Read as double-precision complex then cast
        FeffComplex tmp_dc[nex];
        common::read_pad_complex(fin, npadx, tmp_dc, d.ne);
        for (int i = 0; i < d.ne; ++i) {
            d.phc[i] = std::complex<float>(
                static_cast<float>(tmp_dc[i].real()),
                static_cast<float>(tmp_dc[i].imag()));
        }
    }
    // ck(ne) -- complex
    {
        FeffComplex tmp_dc[nex];
        common::read_pad_complex(fin, npadx, tmp_dc, d.ne);
        for (int i = 0; i < d.ne; ++i) {
            d.ck[i] = std::complex<float>(
                static_cast<float>(tmp_dc[i].real()),
                static_cast<float>(tmp_dc[i].imag()));
        }
    }
    // xk(ne) -- real
    {
        double tmp_d[nex];
        common::read_pad_double(fin, npadx, tmp_d, d.ne);
        for (int i = 0; i < d.ne; ++i) {
            d.xk[i] = static_cast<float>(tmp_d[i]);
        }
    }

    // Allocate per-path arrays
    d.index.resize(npx, 0);
    d.nleg.resize(npx, 0);
    d.deg.resize(npx, 0.0f);
    d.reff.resize(npx, 0.0f);
    d.crit.resize(npx, 0.0f);
    d.ipot.resize(npx, std::vector<int>(legtot, 0));
    d.rat.resize(npx, std::vector<std::array<float, 3>>(legtot, {0.0f, 0.0f, 0.0f}));
    d.beta.resize(npx, std::vector<float>(legtot, 0.0f));
    d.eta.resize(npx, std::vector<float>(legtot, 0.0f));
    d.ri.resize(npx, std::vector<float>(legtot, 0.0f));
    d.achi.resize(npx, std::vector<float>(nex, 0.0f));
    d.phchi.resize(npx, std::vector<float>(nex, 0.0f));

    d.nptot = 0;

    // Read paths
    for (int ipath = 0; ipath < npx; ++ipath) {
        if (!std::getline(fin, str)) break;  // EOF
        pos = str.find_first_not_of(" \t");
        if (pos != std::string::npos) str = str.substr(pos);
        if (str.substr(0, 2) != "##") {
            throw std::runtime_error("rdfbin error: expected ## line in " + fname);
        }

        std::istringstream iss(str.substr(2));
        double tmp_reff = 0.0, tmp_crit = 0.0;
        if (!(iss >> d.index[ipath] >> d.nleg[ipath])) {
            throw std::runtime_error("rdfbin error: bad path header in " + fname);
        }
        float deg_tmp = 0.0f;
        {
            double tmp;
            iss >> tmp;
            deg_tmp = static_cast<float>(tmp);
        }
        d.deg[ipath] = deg_tmp;
        iss >> tmp_reff >> tmp_crit;
        d.reff[ipath] = static_cast<float>(tmp_reff / bohr);
        d.crit[ipath] = static_cast<float>(tmp_crit);

        int nl = d.nleg[ipath];
        for (int j = 0; j < nl; ++j) {
            iss >> d.ipot[ipath][j];
        }

        d.nptot++;

        // Read PAD arrays for this path
        // rat(3*nleg) -- stored flat
        {
            int nvals = 3 * nl;
            std::vector<double> tmp(nvals);
            common::read_pad_double(fin, npadx, tmp.data(), nvals);
            for (int j = 0; j < nl; ++j) {
                d.rat[ipath][j][0] = static_cast<float>(tmp[3 * j]);
                d.rat[ipath][j][1] = static_cast<float>(tmp[3 * j + 1]);
                d.rat[ipath][j][2] = static_cast<float>(tmp[3 * j + 2]);
            }
        }
        // beta(nleg)
        {
            std::vector<double> tmp(nl);
            common::read_pad_double(fin, npadx, tmp.data(), nl);
            for (int j = 0; j < nl; ++j)
                d.beta[ipath][j] = static_cast<float>(tmp[j]);
        }
        // eta(nleg)
        {
            std::vector<double> tmp(nl);
            common::read_pad_double(fin, npadx, tmp.data(), nl);
            for (int j = 0; j < nl; ++j)
                d.eta[ipath][j] = static_cast<float>(tmp[j]);
        }
        // ri(nleg)
        {
            std::vector<double> tmp(nl);
            common::read_pad_double(fin, npadx, tmp.data(), nl);
            for (int j = 0; j < nl; ++j)
                d.ri[ipath][j] = static_cast<float>(tmp[j]);
        }
        // achi(ne)
        {
            std::vector<double> tmp(d.ne);
            common::read_pad_double(fin, npadx, tmp.data(), d.ne);
            for (int ie = 0; ie < d.ne; ++ie)
                d.achi[ipath][ie] = static_cast<float>(tmp[ie]);
            // Zero the rest
            for (int ie = d.ne; ie < nex; ++ie)
                d.achi[ipath][ie] = 0.0f;
        }
        // phchi(ne)
        {
            std::vector<double> tmp(d.ne);
            common::read_pad_double(fin, npadx, tmp.data(), d.ne);
            for (int ie = 0; ie < d.ne; ++ie)
                d.phchi[ipath][ie] = static_cast<float>(tmp[ie]);
            for (int ie = d.ne; ie < nex; ++ie)
                d.phchi[ipath][ie] = 0.0f;
        }
    }

    fin.close();
    return d;
}

} // namespace feff::ff2x
