#include "base/flags.hpp"

#include "absl/strings/str_split.h"
#include "glog/logging.h"

namespace principia {
namespace base {

void Flags::Set(std::string_view const flags) {
  for (std::string_view s : absl::StrSplit(flags, ',')) {
    flags_.insert(absl::StrSplit(s, absl::MaxSplits('=', 1)));
  }
}

bool Flags::IsPresent(std::string_view const flag) {
  return flags_.find(std::string(flag)) != flags_.end();
}

std::string const& Flags::Value(std::string_view const flag) {
  auto const it = flags_.find(std::string(flag));
  CHECK(it != flags_.end());
  return it->second;
}

std::map<std::string, std::string> Flags::flags_;

}  // namespace base
}  // namespace principia
