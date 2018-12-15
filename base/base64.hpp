
#pragma once

#include <cstdint>

#include "base/encoder.hpp"

namespace principia {
namespace base {
namespace internal_base64 {

// This function implements RFC 4648 section 5 (base64url).  The encoded text is
// *not* padded.
template<bool null_terminated>
class Base64Encoder : public Encoder<char, null_terminated> {
 public:
  // Encodes |input| into |output|, which must be large enough to hold the
  // encoded form.
  inline void Encode(Array<std::uint8_t const> input,
                     Array<char> output) override;

  // Same as above but the storage is allocated by the callee.
  inline UniqueArray<char> Encode(Array<std::uint8_t const> input) override;

  // Length of the encoded form, in char16_t.
  inline std::int64_t EncodedLength(Array<std::uint8_t const> input) override;

  // Decodes |input| into |output|, which must be large enough to hold the
  // decoded form.
  inline void Decode(Array<char const> input,
                     Array<std::uint8_t> output) override;

  // Same as above but the storage is allocated by the callee.  The input may or
  // may not be null-terminated.
  inline UniqueArray<std::uint8_t> Decode(Array<char const> input) override;

  // Length of the decoded form, in uint8_t.
  inline std::int64_t DecodedLength(Array<char const> input) override;
};

}  // namespace internal_base64

using internal_base64::Base64Encoder;

}  // namespace base
}  // namespace principia

#include "base/base64_body.hpp"
