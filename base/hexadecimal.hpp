
#pragma once

#include <cstdint>

#include "base/encoder.hpp"

namespace principia {
namespace base {
namespace internal_hexadecimal {

template<bool null_terminated>
class HexadecimalEncoder : public Encoder<char, null_terminated> {
 public:
  // The result is upper-case.  Either |input.data <= &output.data[1]| or
  // |&output.data[input.size << 1] <= input.data| must hold, in particular,
  // |input.data == output.data| is valid.  |output.size| must be at least twice
  // |input.size|.  The range
  // [&output.data[input.size << 1], &output.data[output.size][ is left
  // unmodified.
  void Encode(Array<std::uint8_t const> input,
              Array<char> output) override;

  // Same as above but the storage is allocated by the callee.  If
  // |null_terminated| is true a null byte is appended to the encoded form.
  UniqueArray<char> Encode(Array<std::uint8_t const> input) override;

  // Length of the encoded form, in char.
  std::int64_t EncodedLength(Array<std::uint8_t const> input) override;

  // Invalid digits are read as 0.  If |input.size| is odd, the last character
  // of the input is ignored.  Ignores case.  Either |output.data <=
  // &input.data[1]| or |&input.data[input.size & ~1] <= output.data| must hold,
  // in particular, |input.data == output.data| is valid.  |output.size| must be
  // at least |input.size / 2|.  The range
  // [&output[input.size / 2], &output[output.size][ is left unmodified.
  void Decode(Array<char const> input,
              Array<std::uint8_t> output) override;

  // Same as above but the storage is allocated by the callee.
  UniqueArray<std::uint8_t> Decode(Array<char const> input) override;

  // Length of the decoded form, in uint8_t.
  std::int64_t DecodedLength(Array<char const> input) override;
};

}  // namespace internal_hexadecimal

using internal_hexadecimal::HexadecimalEncoder;

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
