// File I/O for PAD-format data files and headers.
// Converted from: src/COMMON/rdpot.f, src/COMMON/rdxsph.f,
//                 src/COMMON/head.f, src/COMMON/rdhead.f, src/COMMON/rdcmt.f

#include "file_io.hpp"
#include "pad_io.hpp"
#include "logging.hpp"
#include "string_utils.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace feff::common {

// ===================================================================
// read_pot  --  replaces Fortran subroutine rdpot
// ===================================================================
PotData read_pot(const std::string& filename) {
    PotData d;

    std::ifstream in(filename, std::ios::in);
    if (!in.is_open()) {
        par::par_stop("read_pot: cannot open " + filename);
    }

    // ---- First line: 9 integers in format 9(1x,i4)
    {
        std::string line;
        std::getline(in, line);
        // Parse 9 integers from the line
        int ntitle_tmp = 0, nph_tmp = 0, npadx_tmp = 0;
        std::sscanf(line.c_str(), " %d %d %d %d %d %d %d %d %d",
                    &ntitle_tmp, &nph_tmp, &npadx_tmp,
                    &d.nohole, &d.ihole, &d.inters, &d.iafolp,
                    &d.jumprm, &d.iunf);
        d.ntitle = ntitle_tmp;
        d.nph    = nph_tmp;
        d.npadx  = npadx_tmp;
    }

    const int nph = d.nph;
    const int npadx = d.npadx;

    // ---- Title lines
    for (int i = 0; i < d.ntitle; ++i) {
        std::string line;
        std::getline(in, line);
        d.title[i] = ltrim(line);
    }

    // ---- 13 miscellaneous doubles via PAD
    {
        double dum[13]{};
        read_pad_double(in, npadx, dum, 13);
        d.rnrmav = dum[0];
        d.xmu    = dum[1];
        d.vint   = dum[2];
        d.rhoint = dum[3];
        d.emu    = dum[4];
        d.s02    = dum[5];
        d.erelax = dum[6];
        d.wp     = dum[7];
        d.ecv    = dum[8];
        d.rs     = dum[9];
        d.xf     = dum[10];
        d.qtotel = dum[11];
        d.totvol = dum[12];
    }

    // ---- imt  (formatted: 20(1x,i4))
    {
        std::string line;
        std::getline(in, line);
        std::istringstream ss(line);
        for (int i = 0; i <= nph; ++i) ss >> d.imt[i];
    }

    // ---- rmt via PAD
    read_pad_double(in, npadx, d.rmt.data(), nph + 1);

    // ---- inrm
    {
        std::string line;
        std::getline(in, line);
        std::istringstream ss(line);
        for (int i = 0; i <= nph; ++i) ss >> d.inrm[i];
    }

    // ---- iz
    {
        std::string line;
        std::getline(in, line);
        std::istringstream ss(line);
        for (int i = 0; i <= nph; ++i) ss >> d.iz[i];
    }

    // ---- kappa (1:30)
    {
        std::string line;
        std::getline(in, line);
        std::istringstream ss(line);
        for (int i = 0; i < 30; ++i) ss >> d.kappa[i];
    }

    // ---- PAD double arrays (nph+1 sized)
    read_pad_double(in, npadx, d.rnrm.data(),   nph + 1);
    read_pad_double(in, npadx, d.folp.data(),    nph + 1);
    read_pad_double(in, npadx, d.folpx.data(),   nph + 1);
    read_pad_double(in, npadx, d.xnatph.data(),  nph + 1);
    read_pad_double(in, npadx, d.xion.data(),    nph + 1);

    // ---- Dirac spinor components
    read_pad_double(in, npadx, d.dgc0.data(), 251);
    read_pad_double(in, npadx, d.dpc0.data(), 251);

    const int n_orb = 251 * 30 * (nph + 1);
    d.dgc.resize(n_orb);
    d.dpc.resize(n_orb);
    read_pad_double(in, npadx, d.dgc.data(), n_orb);
    read_pad_double(in, npadx, d.dpc.data(), n_orb);

    const int n_ad = 10 * 30 * (nph + 1);
    d.adgc.resize(n_ad);
    d.adpc.resize(n_ad);
    read_pad_double(in, npadx, d.adgc.data(), n_ad);
    read_pad_double(in, npadx, d.adpc.data(), n_ad);

    // ---- Density / potential arrays  251*(nph+1)
    const int n251 = 251 * (nph + 1);
    d.edens.resize(n251);
    d.vclap.resize(n251);
    d.vtot.resize(n251);
    d.edenvl.resize(n251);
    d.vvalgs.resize(n251);
    d.dmag.resize(n251);
    read_pad_double(in, npadx, d.edens.data(),  n251);
    read_pad_double(in, npadx, d.vclap.data(),  n251);
    read_pad_double(in, npadx, d.vtot.data(),   n251);
    read_pad_double(in, npadx, d.edenvl.data(), n251);
    read_pad_double(in, npadx, d.vvalgs.data(), n251);
    read_pad_double(in, npadx, d.dmag.data(),   n251);

    // ---- xnval  30*(nph+1)
    const int n30 = 30 * (nph + 1);
    d.xnval.resize(n30);
    read_pad_double(in, npadx, d.xnval.data(), n30);

    // ---- eorb (1:30)
    read_pad_double(in, npadx, d.eorb.data(), 30);

    // ---- iorb(-4:3, 0:nph)  => 8 values per iph, formatted 8(1x,i2)
    const int n_iorb = 8 * (nph + 1);
    d.iorb.resize(n_iorb);
    for (int iph = 0; iph <= nph; ++iph) {
        std::string line;
        std::getline(in, line);
        std::istringstream ss(line);
        for (int i = 0; i < 8; ++i) {
            ss >> d.iorb[iph * 8 + i];
        }
    }

    // ---- qnrm
    read_pad_double(in, npadx, d.qnrm.data(), nph + 1);

    // ---- xnmues  (lx+1)*(nph+1)
    const int nn = (lx + 1) * (nph + 1);
    d.xnmues.resize(nn);
    read_pad_double(in, npadx, d.xnmues.data(), nn);

    return d;
}

// ===================================================================
// read_xsph  --  replaces Fortran subroutine rdxsph
// ===================================================================
PhaseData read_xsph(const std::string& filename) {
    PhaseData d;

    // Try the supplied path first, fall back to "phase.pad"
    std::ifstream in(filename, std::ios::in);
    if (!in.is_open()) {
        in.open("phase.pad", std::ios::in);
        if (!in.is_open()) {
            par::par_stop("read_xsph: cannot find phase.pad");
        }
    }
    check_file_open(in, filename, "read_xsph");

    // ---- First line: 9 ints + 2 floats
    //      format(9(1x,i4), 2(1x,f10.5))
    {
        std::string line;
        std::getline(in, line);
        d.ixc = 0;
        d.rs  = 0.0;
        d.vint = 0.0;
        // Parse integers then floats
        std::istringstream ss(line);
        ss >> d.nsp >> d.ne >> d.ne1 >> d.ne3 >> d.nph >> d.ihole
           >> d.ik0 >> d.npadx >> d.ixc >> d.rs >> d.vint;
    }

    const int ne    = d.ne;
    const int nph   = d.nph;
    const int npadx = d.npadx;
    const int nsp   = d.nsp;

    // ---- 3 doubles via PAD: rnrmav, xmu, edge
    {
        double dum[3]{};
        read_pad_double(in, npadx, dum, 3);
        d.rnrmav = dum[0];
        d.xmu    = dum[1];
        d.edge   = dum[2];
    }

    // ---- em(ne)  complex
    d.em.resize(ne);
    read_pad_complex(in, npadx, d.em.data(), ne);

    // ---- eref(ne, nsp)  read as flat ne*nsp complex block
    {
        std::vector<FeffComplex> temp(ne * nsp);
        read_pad_complex(in, npadx, temp.data(), ne * nsp);
        d.eref.resize(ne * nsp);
        // Fortran order: inner loop ie, outer loop isp
        int ii = 0;
        for (int isp = 0; isp < nsp; ++isp) {
            for (int ie = 0; ie < ne; ++ie) {
                d.eref[isp * ne + ie] = temp[ii++];
            }
        }
    }

    // ---- Per-potential: lmax0, iz, potlbl, then phase shifts
    // ph stored as ph(ie, ll, isp, iph) in Fortran
    // We flatten to: ph[ iph * (nsp * (2*ltot+1) * ne) +
    //                     isp * ((2*ltot+1) * ne) +
    //                     (ll+ltot) * ne + ie ]
    const int ltot2p1 = 2 * ltot + 1;
    d.ph.resize(static_cast<size_t>(ne) * ltot2p1 * nsp * (nph + 1),
                FeffComplex(0.0, 0.0));
    d.lmax.resize(static_cast<size_t>(ne) * (nph + 1), 0);

    std::vector<int> lmax0(nph + 1);

    for (int iph = 0; iph <= nph; ++iph) {
        // Read lmax0, iz, potlbl  format(2(1x,i3), 1x, a6)
        std::string line;
        std::getline(in, line);
        {
            std::istringstream ss(line);
            ss >> lmax0[iph] >> d.iz[iph];
            // potlbl is 6 chars after the two integers
            size_t pos = 0;
            // Find position after second integer
            // Skip leading whitespace + first int + whitespace + second int
            int count = 0;
            bool in_num = false;
            for (size_t i = 0; i < line.size(); ++i) {
                bool is_ws = (line[i] == ' ' || line[i] == '\t');
                if (!is_ws && !in_num) { in_num = true; ++count; }
                if (is_ws && in_num) { in_num = false; }
                if (count == 2 && is_ws) {
                    // We've passed 2 numbers, skip the space
                    pos = i + 1;
                    break;
                }
            }
            if (pos < line.size()) {
                size_t end = std::min(pos + 6, line.size());
                d.potlbl[iph] = line.substr(pos, end - pos);
            }
        }

        // Read phase shifts for each spin
        for (int isp = 0; isp < nsp; ++isp) {
            const int ii_count = ne * (2 * lmax0[iph] + 1);
            std::vector<FeffComplex> temp(ii_count);
            read_pad_complex(in, npadx, temp.data(), ii_count);

            // Unpack: Fortran loops ie=1,ne then ll=-lmax0..lmax0
            int idx = 0;
            for (int ie = 0; ie < ne; ++ie) {
                for (int ll = -lmax0[iph]; ll <= lmax0[iph]; ++ll) {
                    // ph(ie, ll, isp, iph)
                    size_t flat = static_cast<size_t>(iph) * nsp * ltot2p1 * ne
                                + static_cast<size_t>(isp) * ltot2p1 * ne
                                + static_cast<size_t>(ll + ltot) * ne
                                + ie;
                    d.ph[flat] = temp[idx++];
                }
            }
        }
    }

    // ---- rkk(ne, 8, nsp)  read as flat ne*8*nsp block
    {
        std::vector<FeffComplex> temp(ne * 8 * nsp);
        read_pad_complex(in, npadx, temp.data(), ne * 8 * nsp);
        d.rkk.resize(ne * 8 * nsp);
        // Fortran order: inner ie, then kdif, then isp
        int ii = 0;
        for (int isp = 0; isp < nsp; ++isp) {
            for (int kdif = 0; kdif < 8; ++kdif) {
                for (int ie = 0; ie < ne; ++ie) {
                    d.rkk[isp * 8 * ne + kdif * ne + ie] = temp[ii++];
                }
            }
        }
    }

    // ---- Compute effective lmax per (ie, iph) and lmaxp1
    static constexpr double phmin = 1.0e-7;
    d.lmaxp1 = 0;
    for (int iph = 0; iph <= nph; ++iph) {
        for (int ie = 0; ie < ne; ++ie) {
            int lm = 0;
            for (int il = lmax0[iph]; il >= 0; --il) {
                lm = il;
                // Check spin 0
                size_t idx0 = static_cast<size_t>(iph) * nsp * ltot2p1 * ne
                            + static_cast<size_t>(0) * ltot2p1 * ne
                            + static_cast<size_t>(il + ltot) * ne + ie;
                bool nonzero = std::abs(std::sin(d.ph[idx0])) > phmin;
                if (!nonzero && nsp > 1) {
                    size_t idx1 = static_cast<size_t>(iph) * nsp * ltot2p1 * ne
                                + static_cast<size_t>(nsp - 1) * ltot2p1 * ne
                                + static_cast<size_t>(il + ltot) * ne + ie;
                    nonzero = std::abs(std::sin(d.ph[idx1])) > phmin;
                }
                if (nonzero) break;
            }
            d.lmax[iph * ne + ie] = lm;
            if (lm + 1 > d.lmaxp1) d.lmaxp1 = lm + 1;
        }
    }

    return d;
}

// ===================================================================
// read_header  --  replaces Fortran subroutine rdhead
// ===================================================================
void read_header(std::istream& in, std::vector<std::string>& headers) {
    headers.clear();
    std::string line;
    while (std::getline(in, line)) {
        // End-of-header marker: "--------" starting at position 3 (Fortran col 4:11)
        if (line.size() >= 11 && line.substr(3, 8) == "--------") {
            return;
        }
        headers.push_back(line);
    }
}

// ===================================================================
// write_header  --  replaces Fortran subroutine wthead
// ===================================================================
void write_header(std::ostream& out, const std::vector<std::string>& headers) {
    for (const auto& h : headers) {
        std::string trimmed = rtrim(h);
        out << trimmed << "\n";
    }
}

// ===================================================================
// skip_comments  --  replaces Fortran subroutine rdcmt
// ===================================================================
void skip_comments(std::istream& in, const std::string& comment_chars) {
    // Read lines; if the first non-whitespace character is one of
    // comment_chars, skip it.  Otherwise back up and return.
    std::string line;
    while (in.good()) {
        auto pos = in.tellg();
        if (!std::getline(in, line)) break;

        // Get first non-whitespace character
        bool is_comment = false;
        for (char c : line) {
            if (c == ' ' || c == '\t') continue;
            // Check if this char is in comment_chars
            if (comment_chars.find(c) != std::string::npos) {
                is_comment = true;
            }
            break;
        }

        if (!is_comment) {
            // Back up to before this line
            in.seekg(pos);
            return;
        }
    }
}

// ===================================================================
// make_header  --  replaces Fortran subroutine sthead
// ===================================================================
void make_header(int& ntitle, std::array<std::string, nheadx>& title,
                 int nph,
                 const int* iz, const double* rmt, const double* rnrm,
                 const double* xion, int ihole, int ixc,
                 double vr0, double vi0, double gamach,
                 double xmu, double xf, double vint, double rs,
                 int lreal, double rgrd) {

    static const char* shole[] = {
        "no hole",   "K  shell",  "L1 shell",  "L2 shell",
        "L3 shell",  "M1 shell",  "M2 shell",  "M3 shell",
        "M4 shell",  "M5 shell",  "N1 shell",  "N2 shell",
        "N3 shell",  "N4 shell",  "N5 shell",  "N6 shell",
        "N7 shell",  "O1 shell",  "O2 shell",  "O3 shell",
        "O4 shell",  "O5 shell",  "O6 shell",  "O7 shell",
        "P1 shell",  "P2 shell",  "P3 shell",  "P4 shell",
        "P5 shell",  "R1 shell"
    };
    static const char* sout[] = {
        "H-L exch", "D-H exch", "Gd state", "DH - HL ",
        "DH + HL ", "val=s+d ", "sigmd(r)", "sigmd=c "
    };

    const std::string vfeff_str = "Feff8L (EXAFS) 0.1";

    // First title line: merge user title with version string
    char buf[81]{};
    if (ntitle >= 1 && rtrim(title[0]).size() > 1) {
        std::string t1 = title[0].substr(0, 45);
        std::snprintf(buf, sizeof(buf), "%-45s  %30s", t1.c_str(), vfeff_str.c_str());
    } else {
        std::snprintf(buf, sizeof(buf), "%47s%30s", "", vfeff_str.c_str());
    }
    title[0] = buf;
    int nstor = 1;

    // Remove empty title lines
    for (int it = 1; it < ntitle; ++it) {
        if (rtrim(title[it]).size() > 1) {
            title[nstor] = title[it];
            ++nstor;
        }
    }
    ntitle = nstor;

    // Absorber line
    if (xion[0] != 0.0) {
        std::snprintf(buf, sizeof(buf),
                      "Abs   Z=%2d Rmt=%6.3f Rnm=%6.3f Ion=%5.2f %s",
                      iz[0], rmt[0] * bohr, rnrm[0] * bohr, xion[0],
                      shole[ihole]);
    } else {
        std::snprintf(buf, sizeof(buf),
                      "Abs   Z=%2d Rmt=%6.3f Rnm=%6.3f %s",
                      iz[0], rmt[0] * bohr, rnrm[0] * bohr,
                      shole[ihole]);
    }
    title[ntitle++] = buf;

    // RPHASES / RSIGMA / RGRID line
    if (lreal >= 1 || std::abs(rgrd - 0.05) > 1.0e-5) {
        std::string s1;
        if (lreal > 1)      s1 = "RPHASES";
        else if (lreal == 1) s1 = "RSIGMA";

        std::string s2;
        if (std::abs(rgrd - 0.05) > 1.0e-5) {
            char tmp[20];
            std::snprintf(tmp, sizeof(tmp), "  RGRID%7.4f", rgrd);
            s2 = tmp;
        }
        title[ntitle++] = s1 + s2;
    }

    // Potential lines
    for (int iph = 1; iph <= nph; ++iph) {
        if (xion[iph] != 0.0) {
            std::snprintf(buf, sizeof(buf),
                          "Pot%2d Z=%2d Rmt=%6.3f Rnm=%6.3f Ion=%5.2f",
                          iph, iz[iph], rmt[iph] * bohr, rnrm[iph] * bohr,
                          xion[iph]);
        } else {
            std::snprintf(buf, sizeof(buf),
                          "Pot%2d Z=%2d Rmt=%6.3f Rnm=%6.3f",
                          iph, iz[iph], rmt[iph] * bohr, rnrm[iph] * bohr);
        }
        title[ntitle++] = buf;
    }

    // Gamma / exchange line
    // Fortran: ntitle = ntitle + 1; write(title(ntitle),...) gamach
    if (std::abs(vi0) > 1.0e-8 || std::abs(vr0) > 1.0e-8) {
        std::snprintf(buf, sizeof(buf),
                      "Gam_ch=%9.3E %s Vi=%10.3E Vr=%10.3E",
                      gamach * hart, sout[ixc], vi0 * hart, vr0 * hart);
    } else {
        std::snprintf(buf, sizeof(buf),
                      "Gam_ch=%9.3E %s",
                      gamach * hart, sout[ixc]);
    }
    title[ntitle++] = buf;

    // Mu / kf / Vint / Rs line
    // Fortran: ntitle = ntitle + 1; write(title(ntitle),200) xmu*hart, ...
    std::snprintf(buf, sizeof(buf),
                  "Mu=%10.3E kf=%9.3E Vint=%10.3E Rs_int=%6.3f",
                  xmu * hart, xf / bohr, vint * hart, rs);
    title[ntitle++] = buf;
}

} // namespace feff::common
