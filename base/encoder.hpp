
#pragma once

#include <cstdint>

#include "base/array.hpp"

namespace principia {
namespace base {
namespace internal_encoder {

// Encodes/decodes an array of bytes to/from an array of Char.  If
// |null_terminated| is true a null Char is appended to the encoded form and
// taken into account in the encoded length.
template<typename Char_, bool null_terminated>
class Encoder {
 public:
  using Char = Char_;

  virtual ~Encoder() = default;

  // Encodes |input| into |output|, which must be large enough to hold the
  // encoded form.
  inline virtual void Encode(Array<std::uint8_t const> input,
                             Array<Char> output) = 0;

  // Same as above but the storage is allocated by the callee.
  inline virtual UniqueArray<Char> Encode(Array<std::uint8_t const> input) = 0;

  // Length of the encoded form, in Char.
  inline virtual std::int64_t
  EncodedLength(Array<std::uint8_t const> input) = 0;

  // Decodes |input| into |output|, which must be large enough to hold the
  // decoded form.  The input may or may not be null-terminated.
  inline virtual void Decode(Array<Char const> input,
                             Array<std::uint8_t> output) = 0;

  // Same as above but the storage is allocated by the callee.
  inline virtual UniqueArray<std::uint8_t> Decode(Array<Char const> input) = 0;

  // Length of the decoded form, in uint8_t.
  inline virtual std::int64_t DecodedLength(Array<Char const> input) = 0;
};

}  // namespace internal_encoder

using internal_encoder::Encoder;

}  // namespace base
}  // namespace principia
