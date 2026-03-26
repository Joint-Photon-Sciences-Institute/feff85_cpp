// Packed ASCII Data (PAD) I/O — converted from src/COMMON/padlib.f
//
// Copyright (c) 1997-2001 Matthew Newville, The University of Chicago
// Copyright (c) 1992-1996 Matthew Newville, University of Washington

#include "pad_io.hpp"
#include "logging.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace feff::common {

using namespace pad_constants;

static constexpr double ten = 10.0;
static constexpr double one = 1.0;
static constexpr double zero = 0.0;
static constexpr double base = static_cast<double>(ibase);
static constexpr double pad_huge = std::pow(10.0, ihuge);
static constexpr double pad_tiny = 1.0 / std::pow(10.0, ihuge);

void sclean(std::string& str) {
    for (size_t i = 0; i < str.size(); ++i) {
        int c = static_cast<unsigned char>(str[i]);
        if (c == 0 || (c >= 10 && c <= 15)) {
            // End-of-line: blank remainder
            for (size_t j = i; j < str.size(); ++j) {
                str[j] = ' ';
            }
            return;
        } else if (c <= 31) {
            str[i] = ' ';
        }
    }
}

void pad_encode(double value, int npack, std::string& str) {
    str.assign(npack, ' ');

    double xsave = std::min(pad_huge, std::max(-pad_huge, value));
    int isgn = (xsave > 0.0) ? 1 : 0;

    double xwork = std::abs(xsave);
    int iexp = 0;

    if (xwork < pad_huge && xwork > pad_tiny) {
        iexp = 1 + static_cast<int>(std::log(xwork) / tenlog);
    } else if (xwork >= pad_huge) {
        iexp = ihuge;
        xwork = one;
    } else {
        // xwork <= tiny
        xwork = zero;
    }

    // Force xwork between ~0.1 and ~1
    xwork = xwork / std::pow(ten, iexp);

    // Normalize
    constexpr double tenth = 0.099999999994;
    while (xwork >= one) {
        xwork *= 0.1;
        iexp += 1;
    }
    while (xwork > 0.0 && xwork <= tenth) {
        xwork *= ten;
        iexp -= 1;
    }

    int itmp = static_cast<int>(ibas2 * xwork);
    str[0] = static_cast<char>(iexp + ioff + ibas2);
    str[1] = static_cast<char>(2 * itmp + isgn + ioff);
    xwork = xwork * ibas2 - itmp;

    if (npack > 2) {
        for (int i = 2; i < npack; ++i) {
            itmp = static_cast<int>(base * xwork + 1.0e-9);
            str[i] = static_cast<char>(itmp + ioff);
            xwork = xwork * base - itmp;
        }
    }

    // Rounding
    if (xwork >= 0.5) {
        int ic = itmp + ioff + 1;
        if (ic <= 126) {
            str[npack - 1] = static_cast<char>(ic);
        } else if (npack >= 2) {
            int jc = static_cast<unsigned char>(str[npack - 2]);
            if (jc < 126) {
                str[npack - 2] = static_cast<char>(jc + 1);
                str[npack - 1] = static_cast<char>(37);
            }
        }
    }
}

double pad_decode(const std::string& str, int npack) {
    if (npack <= 2) return zero;

    int iexp = (static_cast<unsigned char>(str[0]) - ioff) - ibas2;
    int char2 = static_cast<unsigned char>(str[1]) - ioff;
    int isgn_val = (char2 % 2) * 2 - 1;   // -1 or +1
    int itmp = char2 / 2;

    double sum = static_cast<double>(itmp) / (base * base);

    // Sum from high index to low for numerical stability (matches Fortran)
    for (int i = npack - 1; i >= 2; --i) {
        sum += static_cast<double>(static_cast<unsigned char>(str[i]) - ioff)
               / std::pow(base, i + 1);
    }

    return 2.0 * isgn_val * ibase * sum * std::pow(ten, iexp);
}

void write_pad_double(std::ostream& out, int npack, const double* arr, int npts) {
    int js = 0;
    std::string str(maxlen, ' ');
    int mxl = maxlen - npack + 1;

    for (int i = 0; i < npts; ++i) {
        std::string encoded;
        pad_encode(arr[i], npack, encoded);
        for (int k = 0; k < npack; ++k) {
            str[js + k] = encoded[k];
        }
        js += npack;

        if (js >= mxl || i == npts - 1) {
            out << cpadr << str.substr(0, js) << '\n';
            js = 0;
            str.assign(maxlen, ' ');
        }
    }
}

void write_pad_complex(std::ostream& out, int npack, const FeffComplex* arr, int npts) {
    int js = 0;
    std::string str(maxlen, ' ');
    int mxl = maxlen - 2 * npack + 1;
    int np = 2 * npack;

    for (int i = 0; i < npts; ++i) {
        std::string re_enc, im_enc;
        pad_encode(arr[i].real(), npack, re_enc);
        pad_encode(arr[i].imag(), npack, im_enc);

        for (int k = 0; k < npack; ++k) {
            str[js + k] = re_enc[k];
            str[js + npack + k] = im_enc[k];
        }
        js += np;

        if (js >= mxl || i == npts - 1) {
            out << cpadc << str.substr(0, js) << '\n';
            js = 0;
            str.assign(maxlen, ' ');
        }
    }
}

/// Internal: read a line, clean it, return useful length (-1 on EOF, -2 on error)
static int iread_line(std::istream& in, std::string& line) {
    if (!std::getline(in, line)) {
        if (in.eof()) return -1;
        return -2;
    }
    sclean(line);
    // Find last non-blank
    size_t last = line.find_last_not_of(" \t\0");
    if (last == std::string::npos) return 0;
    return static_cast<int>(last + 1);
}

int read_pad_double(std::istream& in, int npack, double* arr, int npts) {
    int ipts = 0;
    std::string line;

    while (ipts < npts) {
        int len = iread_line(in, line);
        if (len < 0) break;  // EOF

        // Trim leading whitespace
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        std::string trimmed = line.substr(start);

        char ctest = trimmed[0];
        std::string data = trimmed.substr(1);
        int data_len = static_cast<int>(data.size());
        int ndline = data_len / npack;

        if (ctest != cpadr || ndline <= 0) {
            char buf[256];
            std::snprintf(buf, sizeof(buf), "Read_PAD error: bad data in double array at ipts=%d/%d ctest='%c'(%d) line='%.60s'",
                          ipts, npts, ctest, (int)(unsigned char)ctest, trimmed.c_str());
            throw std::runtime_error(buf);
        }

        for (int i = 0; i < ndline && ipts < npts; ++i) {
            arr[ipts] = pad_decode(data.substr(i * npack, npack), npack);
            ++ipts;
        }
    }

    return ipts;
}

int read_pad_complex(std::istream& in, int npack, FeffComplex* arr, int npts) {
    int ipts = 0;
    int np = 2 * npack;
    std::string line;

    while (ipts < npts) {
        int len = iread_line(in, line);
        if (len < 0) break;  // EOF

        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        std::string trimmed = line.substr(start);

        char ctest = trimmed[0];
        std::string data = trimmed.substr(1);
        int data_len = static_cast<int>(data.size());
        int ndline = data_len / np;

        if (ctest != cpadc || ndline <= 0) {
            throw std::runtime_error("Read_PAD error: bad data in complex array");
        }

        for (int i = 0; i < ndline && ipts < npts; ++i) {
            double re = pad_decode(data.substr(i * np, npack), npack);
            double im = pad_decode(data.substr(i * np + npack, npack), npack);
            arr[ipts] = FeffComplex(re, im);
            ++ipts;
        }
    }

    return ipts;
}

} // namespace feff::common
