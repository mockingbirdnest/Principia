#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace principia {
namespace base {

//TODO(phl):comment, order.
struct UniqueBytes {
  UniqueBytes();  // An object of size 0.
  template<typename T,
           typename = std::enable_if_t<std::is_integral<T>::value, T>>
  explicit UniqueBytes(T const size);
  template<typename T,
           typename = std::enable_if_t<std::is_integral<T>::value, T>>
  UniqueBytes(std::unique_ptr<std::uint8_t[]> data, T const size);

  UniqueBytes(UniqueBytes&& bytes);
  UniqueBytes& operator=(UniqueBytes&& bytes);

  std::unique_ptr<std::uint8_t[]> data;
  std::int64_t size;
};

// A simple container for a pointer and size, similar to a StringPiece.  |data|
// is not owned.
struct Bytes {
  Bytes();  // An object of size 0.
  Bytes(UniqueBytes const& bytes);  // Implicit, no transfer of ownership.
  template<typename T,
           typename = std::enable_if_t<std::is_integral<T>::value, T>>
  Bytes(std::uint8_t* const data, T const size);

  std::uint8_t* data;
  std::int64_t size;
};

// Performs a deep comparison.
bool operator==(Bytes left, Bytes right);
bool operator==(Bytes left, UniqueBytes const& right);
bool operator==(UniqueBytes const& left, Bytes right);
bool operator==(UniqueBytes const& left, UniqueBytes const& right);

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
