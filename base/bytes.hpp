#pragma once

#include <cstdint>

#include "base/not_null.hpp"

namespace principia {
namespace base {

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.
struct Bytes {
  Bytes();  // An object of size 0.
  Bytes(not_null<std::uint8_t*> const data, int const size);
  not_null<std::uint8_t*> data;
  int size;
};

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
