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

  //TODO(phl):Comment
  void CheckNotNull() const;

  static Bytes New(std::int64_t const size);
  static Bytes New(std::string s);
  static void Delete(Bytes const bytes);

  std::uint8_t* data;
  std::int64_t size;
};

// Performs a deep comparison.
static bool operator==(Bytes left, Bytes right);

// A helper class for UniqueBytes.
struct BytesDeleter {
  void operator()(Bytes* const bytes) const {
  Bytes::Delete(*bytes);
  delete bytes;
    }
};

// An RAII object which takes ownership of the |data|.
using UniqueBytes = std::unique_ptr<Bytes, BytesDeleter>;

}  // namespace base
}  // namespace principia

#include "base/bytes_body.hpp"
