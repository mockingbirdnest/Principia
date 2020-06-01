#include "base/flags.hpp"

#include <map>
#include <set>
#include <string>

#include "glog/logging.h"

namespace principia {
namespace base {

void Flags::Clear() {
  flags().clear();
}

void Flags::Set(std::string_view const name, std::string_view const value) {
  LOG(INFO) << "Setting flag " << name << " = " << value;
  flags().emplace(std::string(name), std::string(value));
}

bool Flags::IsPresent(std::string_view const name) {
  return flags().find(std::string(name)) != flags().end();
}

bool Flags::IsPresent(std::string_view const name,
                      std::string_view const value) {
  auto const pair = flags().equal_range(std::string(name));
  std::set<std::string> values;
  for (auto it = pair.first; it != pair.second; ++it) {
    if (it->second == value) {
      return true;
    }
  }
  return false;
}

std::set<std::string> Flags::Values(std::string_view const name) {
  auto const pair = flags().equal_range(std::string(name));
  std::set<std::string> values;
  for (auto it = pair.first; it != pair.second; ++it) {
    values.insert(it->second);
  }
  return values;
}

std::multimap<std::string, std::string>& Flags::flags() {
  static auto* const flags = new std::multimap<std::string, std::string>();
  return *flags;
}

}  // namespace base
}  // namespace principia
