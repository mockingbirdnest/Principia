export module principia.base.base64;

import std;

import "absl/strings/escaping.h";

import principia.base.array;
import principia.base.encoder;
import principia.base.base64;

namespace principia {
namespace base {
namespace _base64 {
namespace internal {

using namespace principia::base::_array;
using namespace principia::base::_encoder;

export {

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

}  // export
}  // namespace internal

export using internal::Base64Encoder;

}  // namespace _base64
}  // namespace base
}  // namespace principia

#include "base/base64_body.hpp"
