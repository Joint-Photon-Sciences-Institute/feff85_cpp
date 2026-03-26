// String utilities — converted from src/COMMON/str.f, isnum.f, str2dp.f

#include "string_utils.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace feff::common {

std::string ltrim(const std::string& s) {
    auto it = std::find_if(s.begin(), s.end(), [](unsigned char c) {
        return c != ' ' && c != '\t';
    });
    return std::string(it, s.end());
}

std::string rtrim(const std::string& s) {
    auto it = std::find_if(s.rbegin(), s.rend(), [](unsigned char c) {
        return c != ' ' && c != '\t' && c != '\0';
    });
    return std::string(s.begin(), it.base());
}

std::string trim(const std::string& s) {
    return ltrim(rtrim(s));
}

std::string to_upper(const std::string& s) {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return result;
}

std::string to_lower(const std::string& s) {
    std::string result = s;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return result;
}

std::vector<std::string> split_words(const std::string& s, int max_words) {
    // Matches Fortran bwords(): words separated by whitespace/tabs or commas.
    // Consecutive commas with only whitespace between them produce a blank word.
    std::vector<std::string> words;
    bool between = true;
    bool comma_found = true;
    size_t begin_pos = 0;

    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == ' ' || c == '\t') {
            if (!between) {
                words.push_back(s.substr(begin_pos, i - begin_pos));
                between = true;
                comma_found = false;
            }
        } else if (c == ',') {
            if (!between) {
                words.push_back(s.substr(begin_pos, i - begin_pos));
                between = true;
            } else if (comma_found) {
                words.emplace_back(" ");
            }
            comma_found = true;
        } else {
            if (between) {
                between = false;
                begin_pos = i;
            }
        }
        if (max_words > 0 && static_cast<int>(words.size()) >= max_words) {
            return words;
        }
    }

    if (!between) {
        words.push_back(s.substr(begin_pos));
    }

    return words;
}

void replace_tabs(std::string& s) {
    std::replace(s.begin(), s.end(), '\t', ' ');
}

bool is_comment(const std::string& line) {
    if (line.empty() || trim(line).empty()) return true;
    // Find first non-blank character (matches Fortran iscomm behavior)
    std::string trimmed = ltrim(line);
    if (trimmed.empty()) return true;
    char first = trimmed[0];
    return first == ';' || first == '*' || first == '%' || first == '#';
}

bool is_numeric(const std::string& s) {
    // Matches Fortran isnum(): allows digits, sign, decimal, whitespace,
    // and at most one exponent marker (d/D/e/E) and one decimal point.
    const std::string allowed = "deDE.,+- 1234567890";
    int exp_count = 0;
    int dec_count = 0;

    std::string trimmed = trim(s);
    if (trimmed.empty()) return false;

    for (char c : trimmed) {
        size_t pos = allowed.find(c);
        if (pos == std::string::npos) return false;
        // positions 0-3 are d,e,D,E
        if (pos <= 3) exp_count++;
        // position 4 is '.'
        if (pos == 4) dec_count++;
    }

    return exp_count <= 1 && dec_count <= 1;
}

bool parse_double(const std::string& s, double& value) {
    std::string trimmed = trim(s);
    if (trimmed.empty()) return false;

    // Fortran allows 'd' as exponent marker
    std::string normalized = trimmed;
    for (auto& c : normalized) {
        if (c == 'd' || c == 'D') c = 'e';
    }

    try {
        size_t pos;
        value = std::stod(normalized, &pos);
        return pos == normalized.size();
    } catch (...) {
        return false;
    }
}

bool parse_float(const std::string& s, float& value) {
    double dval;
    if (!parse_double(s, dval)) return false;
    value = static_cast<float>(dval);
    return true;
}

bool parse_int(const std::string& s, int& value) {
    std::string trimmed = trim(s);
    if (trimmed.empty()) return false;
    try {
        size_t pos;
        value = std::stoi(trimmed, &pos);
        return pos == trimmed.size();
    } catch (...) {
        return false;
    }
}

} // namespace feff::common
