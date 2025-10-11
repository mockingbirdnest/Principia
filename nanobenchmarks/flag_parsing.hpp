#pragma once

#include <string>
#include <string_view>
#include <vector>

// Adding support for flag types only works using ADL (or by being in
// marshalling.h), so we do this, which is UB.
namespace std {

extern bool AbslParseFlag(std::string_view text,
                          std::vector<double>* flag,
                          std::string* error);

extern std::string AbslUnparseFlag(std::vector<double> const& flag);

}  // namespace std
