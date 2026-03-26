#pragma once
// Packed ASCII Data (PAD) format I/O.
// Converted from: src/COMMON/padlib.f + padlib.h
//
// PAD encodes floating-point numbers as printable ASCII strings for
// portable data files (pot.pad, phase.pad). Each number uses npack
// characters with base-90 encoding.
//
// Copyright (c) 1997-2001 Matthew Newville, The University of Chicago
// Copyright (c) 1992-1996 Matthew Newville, University of Washington

#include <feff/types.hpp>
#include <string>
#include <iostream>

namespace feff::common {

/// PAD format constants (from padlib.h)
namespace pad_constants {
    inline constexpr char cpadr = '!';   // real data line marker
    inline constexpr char cpadc = '$';   // complex data line marker
    inline constexpr char cpadi = '%';   // integer data line marker
    inline constexpr int ibase = 90;     // encoding base
    inline constexpr int ioff = 37;      // ASCII offset
    inline constexpr int ihuge = 38;     // max exponent
    inline constexpr int ibas2 = ibase / 2;  // half-base
    inline constexpr int maxlen = 82;    // max line length
    inline constexpr double tenlog = 2.302585092994045684;  // ln(10)
}

/// Encode a double-precision value to a PAD string of length npack.
void pad_encode(double value, int npack, std::string& str);

/// Decode a PAD string of length npack to a double-precision value.
double pad_decode(const std::string& str, int npack);

/// Write a double-precision array in PAD format.
void write_pad_double(std::ostream& out, int npack, const double* arr, int npts);

/// Write a complex (double) array in PAD format.
void write_pad_complex(std::ostream& out, int npack, const FeffComplex* arr, int npts);

/// Read a double-precision array from PAD format.
/// Returns the number of values actually read.
int read_pad_double(std::istream& in, int npack, double* arr, int npts);

/// Read a complex (double) array from PAD format.
/// Returns the number of values actually read.
int read_pad_complex(std::istream& in, int npack, FeffComplex* arr, int npts);

/// Clean a string of non-printable characters (replaces Fortran sclean).
/// Characters char(0) and char(10)-char(15) terminate the line;
/// other control chars (< 32) are replaced with spaces.
void sclean(std::string& str);

} // namespace feff::common
