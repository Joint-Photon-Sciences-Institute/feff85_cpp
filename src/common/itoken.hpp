#pragma once
// Token identifier for FEFF input file keywords.
// Converted from: src/COMMON/itoken.f

#include <string>

namespace feff::common {

/// Return the integer token ID for a keyword in the given input file context.
/// Returns 0 if the word is not recognized, -1 for "END".
/// The word is matched case-insensitively using the first 4 characters.
int itoken(const std::string& word, const std::string& filename);

} // namespace feff::common
