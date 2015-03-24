#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.
struct Bytes {
  Bytes();  // An object of size 0.
  Bytes(std::uint8_t* const data, std::int64_t const size);

  std::uint8_t* data;
  std::int64_t size;
};

struct UniqueBytes {
  UniqueBytes();  // An object of size 0.
  UniqueBytes(std::int64_t const size);
  UniqueBytes(std::unique_ptr<std::uint8_t[]> data,
              std::int64_t const size);

  std::unique_ptr<std::uint8_t[]> data;
  std::int64_t size;
};

// Performs a deep comparison.
static bool operator==(Bytes left, Bytes right);
static bool operator==(UniqueBytes left, UniqueBytes right);

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
