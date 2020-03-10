#pragma once

#include <map>
#include <set>
#include <string_view>

namespace principia {
namespace base {

class Flags {
 public:
  // Remove all the flags.
  static void Clear();

  // Sets a flag with the given |name| and |value|.
  static void Set(std::string_view name, std::string_view value);

  // Returns true if a flag with the given |name| was set (with any value).
  static bool IsPresent(std::string_view name);

  // Returns true if a flag with the given |name| and |value| was set.
  static bool IsPresent(std::string_view name, std::string_view value);

  // Returns the values associated with the flag with a given |name|.  Returns
  // an empty set if there is no flag with the |name|.
  static std::set<std::string> Values(std::string_view name);

 private:
  static std::multimap<std::string, std::string> flags_;
};

}  // namespace base
}  // namespace principia
