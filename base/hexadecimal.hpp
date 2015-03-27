#pragma once

#include <cstdint>

#include "base/array.hpp"

namespace principia {
namespace base {

// The result is upper-case.  Either |input.data <= &output.data[1]| or
// |&output.data[input.size << 1] <= input.data| must hold, in particular,
// |input.data == output.data| is valid.  |output.size| must be at least twice
// |input.size|.  The range
// [&output.data[input.size << 1], &output.data[output.size][ is left
// unmodified.
void HexadecimalEncode(Array<std::uint8_t const> input,
                       Array<std::uint8_t> output);

// Invalid digits are read as 0.  If |input.size| is odd, the last character of
// the input is ignored.  Ignores case.  Either |output.data <= &input.data[1]|
// or |&input.data[input.size & ~1] <= output.data| must hold, in particular,
// |input.data == output.data| is valid.  |output.size| must be at least
// |input.size / 2|.  The range
// [&output[input.size / 2], &output[output.size][ is left unmodified.
void HexadecimalDecode(Array<std::uint8_t const> input,
                       Array<std::uint8_t> output);

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
