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
  template<typename T>
  Bytes(std::uint8_t* const data, T const size);

  std::uint8_t* data;
  std::int64_t size;
};

//TODO(phl):comment.
struct UniqueBytes {
  UniqueBytes();  // An object of size 0.
  template<typename T>
  UniqueBytes(T const size);
  template<typename T>
  UniqueBytes(std::unique_ptr<std::uint8_t[]> data,
              T const size);
  ~UniqueBytes();

  std::uint8_t* data;
  std::int64_t size;
};

// Performs a deep comparison.
bool operator==(Bytes left, Bytes right);
bool operator==(Bytes left, UniqueBytes right);
bool operator==(UniqueBytes left, Bytes right);
bool operator==(UniqueBytes left, UniqueBytes right);

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
