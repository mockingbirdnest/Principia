#include "base/flags.hpp"

#include <functional>
#include <map>
#include <set>
#include <string>
#include <string_view>

#include "absl/log/check.h"
#include "absl/log/log.h"

namespace principia {
namespace base {
namespace _flags {
namespace internal {

void Flags::Clear() {
  flags().clear();
}

void Flags::Set(std::string_view const name, std::string_view const value) {
  LOG(INFO) << "Setting flag " << name << " = " << value;
  flags().emplace(name, value);
}

bool Flags::IsPresent(std::string_view const name) {
  return flags().contains(name);
}

bool Flags::IsPresent(std::string_view const name,
                      std::string_view const value) {
  auto const pair = flags().equal_range(name);
  for (auto it = pair.first; it != pair.second; ++it) {
    if (it->second == value) {
      return true;
    }
  }
  return false;
}

std::set<std::string> Flags::Values(std::string_view const name) {
  auto const pair = flags().equal_range(name);
  std::set<std::string> values;
  for (auto it = pair.first; it != pair.second; ++it) {
    values.insert(it->second);
  }
  return values;
}

std::multimap<std::string, std::string, std::less<>>& Flags::flags() {
  static auto* const flags =
      new std::multimap<std::string, std::string, std::less<>>();
  return *flags;
}

}  // namespace internal
}  // namespace _flags
}  // namespace base
}  // namespace principia
