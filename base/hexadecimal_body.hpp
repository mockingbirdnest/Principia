
#pragma once

#include <cstdint>
#include <cstring>

#include "base/hexadecimal.hpp"
#include "glog/logging.h"

namespace principia {
namespace base {

static char const byte_to_hexadecimal_digits[] =
    "000102030405060708090A0B0C0D0E0F101112131415161718191A1B1C1D1E1F2021222324"
    "25262728292A2B2C2D2E2F303132333435363738393A3B3C3D3E3F40414243444546474849"
    "4A4B4C4D4E4F505152535455565758595A5B5C5D5E5F606162636465666768696A6B6C6D6E"
    "6F707172737475767778797A7B7C7D7E7F808182838485868788898A8B8C8D8E8F90919293"
    "9495969798999A9B9C9D9E9FA0A1A2A3A4A5A6A7A8A9AAABACADAEAFB0B1B2B3B4B5B6B7B8"
    "B9BABBBCBDBEBFC0C1C2C3C4C5C6C7C8C9CACBCCCDCECFD0D1D2D3D4D5D6D7D8D9DADBDCDD"
    "DEDFE0E1E2E3E4E5E6E7E8E9EAEBECEDEEEFF0F1F2F3F4F5F6F7F8F9FAFBFCFDFEFF";

#if defined(SKIP_48) || defined(SKIP_26) || defined(SKIP_7)
#error SKIP_* already defined
#else
#define SKIP_7 0, 0, 0, 0, 0, 0, 0
#define SKIP_26 SKIP_7, SKIP_7, SKIP_7, 0, 0, 0, 0, 0
#define SKIP_48 SKIP_26, SKIP_7, SKIP_7, SKIP_7, 0

static std::uint8_t const hexadecimal_digits_to_nibble[256] = {
    SKIP_48, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    SKIP_7, '\xa', '\xb', '\xc', '\xd', '\xe', '\xf',
    SKIP_26, '\xa', '\xb', '\xc', '\xd', '\xe', '\xf'};

#undef SKIP_7
#undef SKIP_26
#undef SKIP_48
#endif

void HexadecimalEncode(Array<std::uint8_t const> input,
                       Array<std::uint8_t> output) {
  CHECK_NOTNULL(input.data);
  CHECK_NOTNULL(output.data);
  // We iterate backward.
  // |input <= &output[1]| is still valid because we write two bytes of output
  // from reading one byte of input, so output[1] and output[0] are written
  // after reading input[0].  Greater values of |output| would
  // overwrite input data before it is read, unless there is no overlap, i.e.,
  // |&output[input_size << 1] <= input|.
  CHECK(input.data <= &output.data[1] ||
        &output.data[input.size << 1] <= input.data) << "bad overlap";
  CHECK_GE(output.size, input.size << 1) << "output too small";
  // We want the result to start at |output.data[0]|.
  output.data = output.data + ((input.size - 1) << 1);
  input.data = input.data + input.size - 1;
  for (std::uint8_t const* const input_rend = input.data - input.size;
       input.data != input_rend;
       --input.data, output.data -= 2) {
    std::memcpy(output.data, &byte_to_hexadecimal_digits[*input.data << 1], 2);
  }
}
UniqueArray<std::uint8_t> HexadecimalEncode(Array<std::uint8_t const> input,
                                            bool const null_terminated) {
  base::UniqueArray<std::uint8_t> output((input.size << 1) +
                                         (null_terminated ? 1 : 0));
  if (output.size > 0) {
    base::HexadecimalEncode(input, output.get());
  }
  if (null_terminated) {
    output.data[output.size - 1] = 0;
  }
  return output;
}

void HexadecimalDecode(Array<std::uint8_t const> input,
                       Array<std::uint8_t> output) {
  CHECK_NOTNULL(input.data);
  CHECK_NOTNULL(output.data);
  input.size &= ~1;
  // |output <= &input[1]| is still valid because we write one byte of output
  // from reading two bytes of input, so output[0] is written after reading
  // input[0] and input[1].  Greater values of |output| would overwrite input
  // data before it is read, unless there is no overlap, i.e.,
  // |&input[input_size] <= output|.
  CHECK(output.data <= &input.data[1] ||
        &input.data[input.size] <= output.data) << "bad overlap";
  CHECK_GE(output.size, input.size / 2) << "output too small";
  for (std::uint8_t const* const input_end = input.data + input.size;
       input.data != input_end;
       input.data += 2, ++output.data) {
    *output.data = (hexadecimal_digits_to_nibble[*input.data] << 4) |
                   hexadecimal_digits_to_nibble[*(input.data + 1)];
  }
}

UniqueArray<std::uint8_t> HexadecimalDecode(Array<std::uint8_t const> input) {
  UniqueArray<std::uint8_t> output(input.size >> 1);
  if (output.size > 0) {
    HexadecimalDecode({ input.data, input.size & ~1 }, output.get());
  }
  return output;
}

}  // namespace base
}  // namespace principia
