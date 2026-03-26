#pragma once
// String utilities for parsing and text manipulation.
// Converted from: src/COMMON/str.f, src/COMMON/isnum.f, src/COMMON/str2dp.f

#include <string>
#include <vector>

namespace feff::common {

/// Remove leading whitespace (tabs and spaces). Replaces Fortran triml().
std::string ltrim(const std::string& s);

/// Remove trailing whitespace. Replaces checking via istrln().
std::string rtrim(const std::string& s);

/// Remove leading and trailing whitespace.
std::string trim(const std::string& s);

/// Convert string to uppercase. Replaces Fortran upper().
std::string to_upper(const std::string& s);

/// Convert string to lowercase. Replaces Fortran lower().
std::string to_lower(const std::string& s);

/// Split string into words separated by whitespace or commas.
/// Replaces Fortran bwords(). max_words=0 means no limit.
std::vector<std::string> split_words(const std::string& s, int max_words = 0);

/// Replace all tab characters with spaces. Replaces Fortran untab().
void replace_tabs(std::string& s);

/// Test if line is a comment or blank (starts with ;*%#, or is empty).
/// Replaces Fortran iscomm().
bool is_comment(const std::string& line);

/// Test if string can represent a number. Replaces Fortran isnum().
/// Allows digits, sign, decimal point, and at most one 'd'/'e' exponent.
bool is_numeric(const std::string& s);

/// Parse string to double. Replaces Fortran str2dp().
/// Returns true on success, sets value. Returns false on parse error.
bool parse_double(const std::string& s, double& value);

/// Parse string to float. Replaces Fortran str2re().
bool parse_float(const std::string& s, float& value);

/// Parse string to int. Replaces Fortran str2in().
bool parse_int(const std::string& s, int& value);

} // namespace feff::common
