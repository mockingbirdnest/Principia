#include "nanobenchmarks/flags.hpp"

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"

namespace std {

bool AbslParseFlag(std::string_view const text,
                   std::vector<double>* const flag,
                   std::string* const error) {
  flag->clear();
  for (absl::string_view const element : absl::StrSplit(text, ',')) {
    if (!absl::ParseFlag(element, &flag->emplace_back(), error)) {
      return false;
    }
  }
  return true;
}

std::string AbslUnparseFlag(std::vector<double> const& flag) {
  return absl::StrJoin(flag, ",");
}

}  // namespace std

