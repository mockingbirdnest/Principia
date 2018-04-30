
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
inline void HexadecimalEncode(Array<std::uint8_t const> input,
                              Array<std::uint8_t> output);

// Same as above but the storage is allocated by the callee.  If
// |null_terminated| is true a null byte is appended to the encoded form.
inline UniqueArray<std::uint8_t> HexadecimalEncode(
    Array<std::uint8_t const> input,
    bool null_terminated);

// Invalid digits are read as 0.  If |input.size| is odd, the last character of
// the input is ignored.  Ignores case.  Either |output.data <= &input.data[1]|
// or |&input.data[input.size & ~1] <= output.data| must hold, in particular,
// |input.data == output.data| is valid.  |output.size| must be at least
// |input.size / 2|.  The range
// [&output[input.size / 2], &output[output.size][ is left unmodified.
inline void HexadecimalDecode(Array<std::uint8_t const> input,
                              Array<std::uint8_t> output);

// Same as above but the storage is allocated by the callee.  The input may or
// may not be null-terminated.
inline UniqueArray<std::uint8_t> HexadecimalDecode(
    Array<std::uint8_t const> input);

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
