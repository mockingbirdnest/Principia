#pragma once

#include <cstdint>

#include "base/encoder.hpp"

namespace principia {
namespace base {
namespace _base64 {
namespace internal {

using namespace principia::base::_array;

// This function implements RFC 4648 section 5 (base64url).  The encoded text is
// *not* padded.
template<bool null_terminated>
class Base64Encoder : public Encoder<char, null_terminated> {
 public:
  void Encode(Array<std::uint8_t const> input,
              Array<char> output) override;

  UniqueArray<char> Encode(Array<std::uint8_t const> input) override;

  std::int64_t EncodedLength(Array<std::uint8_t const> input) override;

  void Decode(Array<char const> input,
              Array<std::uint8_t> output) override;

  UniqueArray<std::uint8_t> Decode(Array<char const> input) override;

  std::int64_t DecodedLength(Array<char const> input) override;
};

}  // namespace internal

using internal::Base64Encoder;

}  // namespace _base64
}  // namespace base
}  // namespace principia

namespace principia::base {
using namespace principia::base::_base64;
}  // namespace principia::base

#include "base/base64_body.hpp"
