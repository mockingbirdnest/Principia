#pragma once

#include <cstdint>
#include <string>

namespace principia {
namespace base {

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.
struct Bytes {
  Bytes();  // An object of size 0.
  Bytes(std::uint8_t* const data, std::int64_t const size);

  //TODO(phl):Comment
  void CheckNotNull() const;

  static Bytes New(std::int64_t const size);
  static Bytes New(std::string s);
  static void Delete(Bytes const bytes);

  std::uint8_t* data;
  std::int64_t size;
};

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
