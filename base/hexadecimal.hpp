#pragma once

#include <stddef.h>
#include <stdint.h>

namespace principia {
namespace base {

// The result is upper-case.  Either |&input[input_size] <= output| or
// |&output[input_size] <= input| must hold.  |output_size| must be at least
// twice |input_size|.
void HexadecimalEncode(uint8_t const* input, size_t input_size,
                       uint8_t* output, size_t output_size);

// Invalid digits are read as 0.  If |input_size| is odd, the last
// character of the input is ignored.  Ignores case.
// Either |output <= &input[1]| or |&input[input_size & ~1] <= output| must
// hold, in particular, |input == output| is valid.  |output_size| must be at
// least |input_size / 2|.
void HexadecimalDecode(uint8_t const* input, size_t input_size,
                       uint8_t* output, size_t output_size);

}  // namespace base
}  // namespace principia

#include "base/hexadecimal_body.hpp"
