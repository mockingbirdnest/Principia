// See macos_allocator_replacement for an explanation of the pattern used in
// this file.

#pragma once

#include <iomanip>
#include <string>

#include "base/macros.hpp"  // 🧙 For OS_MACOSX.

#if !OS_MACOSX
#error Only include this file for macOS
#endif

namespace principia {
namespace std {

using namespace ::std;

// Apple has an oldish libc++.  This implementation was copied from
// https://github.com/Norgg/eggsperimental_filesystem.
namespace filesystem {

class path {
 public:
  path() = default;
  path(std::string const& source) : native_(source) {};  // NOLINT
  path(char const* source) : native_(source) {};         // NOLINT
  path(path const&) = default;
  ~path() = default;

  path& operator/=(path const& source) {
    native_ += '/';
    // Nah, we're not stripping redundant slashes or anything.
    native_ += source;
    return *this;
  }

  path& operator+=(path const& source) {
    native_ += source;
    return *this;
  }

  std::string const& native() const {
    return native_;
  }

  operator std::string() const {
    return native_;
  }

  std::string string() const {
    return native_;
  }

  // Doesn't actually replace anything.
  path& replace_extension(path const& replacement) {
    native_ += '.';
    native_ += replacement;
    return *this;
  }

  path stem() const {
    auto const pos1 = native_.find_last_of('/');
    auto const pos2 = native_.find_last_of('.');
    return path(native_.substr(pos1 + 1, pos2 - pos1 - 1));
  }

  path extension() const {
    auto const pos = native_.find_last_of('.');
    return path(native_.substr(pos + 1));
  }

  path parent_path() const {
    auto const pos = native_.find_last_of('/');
    return path(native_.substr(0, pos));
  }

 private:
  std::string native_;
};

inline path operator/(path const& left, path const& right) {
  return path(left) /= right;
}

template<typename Character, typename Traits>
std::basic_ostream<Character, Traits>& operator<<(
    std::basic_ostream<Character, Traits>& os,
    const path& p) {
  return os << std::quoted((string)p);
}

}  // namespace filesystem
}  // namespace std
}  // namespace principia
