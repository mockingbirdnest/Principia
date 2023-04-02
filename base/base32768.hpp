#pragma once

#include <cstdint>

#include "base/encoder.hpp"

namespace principia {
namespace base {
namespace _base32768 {
namespace internal {

using namespace principia::base::_array;
using namespace principia::base::_encoder;

// An encoder the base 32768 encoding defined by
// https://github.com/qntm/base32768.  This is a complete reimplementation in
// C++.
template<bool null_terminated>
class Base32768Encoder : public Encoder<char16_t, null_terminated> {
 public:
  void Encode(Array<std::uint8_t const> input,
              Array<char16_t> output) override;

  UniqueArray<char16_t> Encode(Array<std::uint8_t const> input) override;

  std::int64_t EncodedLength(Array<std::uint8_t const> input) override;

  void Decode(Array<char16_t const> input,
              Array<std::uint8_t> output) override;

  UniqueArray<std::uint8_t> Decode(Array<char16_t const> input) override;

  std::int64_t DecodedLength(Array<char16_t const> input) override;
};

}  // namespace internal

using internal::Base32768Encoder;

}  // namespace _base32768
}  // namespace base
}  // namespace principia

#include "base/base32768_body.hpp"
