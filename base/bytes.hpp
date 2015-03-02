#pragma once

#include <cstdint>

#include "base/not_null.hpp"

namespace principia {
namespace base {

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.
struct Bytes {
  Bytes(not_null<std::uint8_t const*> const data, int const size);
  not_null<std::uint8_t const*> const data;
  int const size;

  // A |Bytes| object of size 0.
  static Bytes const Null;

 private:
  // Placeholder for an empty |Bytes| object.
  static std::uint8_t null_data_;
};

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
