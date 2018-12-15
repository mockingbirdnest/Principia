
#pragma once

#include <cstdint>

#include "base/encoder.hpp"

// This file implements the base 32768 encoding defined by
// https://github.com/qntm/base32768.  This is a complete reimplementation in
// C++.

namespace principia {
namespace base {
namespace internal_base32768 {

template<bool null_terminated>
class Base32768Encoder : public Encoder<char16_t, null_terminated> {
 public:
  // Encodes |input| into |output|, which must be large enough to hold the
  // encoded form.
  void Encode(Array<std::uint8_t const> input,
              Array<char16_t> output) override;

  // Same as above but the storage is allocated by the callee.  If
  // |null_terminated| is true a null character is appended to the encoded form.
  UniqueArray<char16_t> Encode(Array<std::uint8_t const> input) override;

  // Length of the encoded form, in char16_t.
  std::int64_t EncodedLength(Array<std::uint8_t const> input) override;

  // Decodes |input| into |output|, which must be large enough to hold the
  // decoded form.
  void Decode(Array<char16_t const> input,
              Array<std::uint8_t> output) override;

  // Same as above but the storage is allocated by the callee.  The input may or
  // may not be null-terminated.
  UniqueArray<std::uint8_t> Decode(Array<char16_t const> input) override;

  // Length of the decoded form, in uint8_t.
  std::int64_t DecodedLength(Array<char16_t const> input) override;
};

}  // namespace internal_base32768

using internal_base32768::Base32768Encoder;

}  // namespace base
}  // namespace principia

#include "base/base32768_body.hpp"
