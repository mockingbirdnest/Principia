#pragma once

#include <map>
#include <string_view>

namespace principia {
namespace base {

class Flags {
 public:
  // |flags| must be a comma-separated string of elements that may contain an
  // equal sign thus:
  //   a=b,c=,d
  static void Set(std::string_view flags);

  // Returns true if the given |flag| was specified as a key in the string
  // passed to Set.
  static bool IsPresent(std::string_view flag);

  // Returns the value associated with the key |flag|, or an empty string if no
  // value was specified.
  static std::string const& Value(std::string_view flag);

 private:
  static std::map<std::string, std::string> flags_;
};

}  // namespace base
}  // namespace principia
