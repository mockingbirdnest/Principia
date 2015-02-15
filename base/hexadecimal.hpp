#pragma once

#include <stdint.h>

namespace principia {
namespace base {

// The result is upper-case.  Either |input <= &output[1]| or
// |&output[input_size << 1] <= input| must hold, in particular,
// |input == output| is valid.  |output_size| must be at least twice
// |input_size|.  The range [&output[input_size << 1], &output[output_size])
// is left unmodified.
void HexadecimalEncode(uint8_t const* input, int64_t const input_size,
                       uint8_t* output, int64_t const output_size);

// Invalid digits are read as 0.  If |input_size| is odd, the last
// character of the input is ignored.  Ignores case.
// Either |output <= &input[1]| or |&input[input_size & ~1] <= output| must
// hold, in particular, |input == output| is valid.  |output_size| must be at
// least |input_size / 2|.  The range
// [&output[input_size / 2], &output[output_size]) is left unmodified.
void HexadecimalDecode(uint8_t const* input, int64_t input_size,
                       uint8_t* output, int64_t const output_size);

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
