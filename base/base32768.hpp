
#pragma once

#include <cstdint>

#include "base/encoder.hpp"

namespace principia {
namespace base {
namespace internal_base32768 {

// An encoder the base 32768 encoding defined by
// https://github.com/qntm/base32768.  This is a complete reimplementation in
// C++.
template<bool null_terminated>
class Base32768Encoder : public Encoder<char16_t, null_terminated> {
 public:
  inline void Encode(Array<std::uint8_t const> input,
                     Array<char16_t> output) override;

  inline UniqueArray<char16_t> Encode(Array<std::uint8_t const> input) override;

  inline std::int64_t EncodedLength(Array<std::uint8_t const> input) override;

  inline void Decode(Array<char16_t const> input,
                     Array<std::uint8_t> output) override;

  inline UniqueArray<std::uint8_t> Decode(Array<char16_t const> input) override;

  inline std::int64_t DecodedLength(Array<char16_t const> input) override;
};

}  // namespace internal_base32768

using internal_base32768::Base32768Encoder;

}  // namespace base
}  // namespace principia

#include "base/base32768_body.hpp"
