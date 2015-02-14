#pragma once

#include <string.h>

#include "base/hexadecimal.hpp"

namespace principia {
namespace base {

static char const kByteToHexadecimalDigits[] = 
"000102030405060708090A0B0C0D0E0F101112131415161718191A1B1C1D1E1F202122232425262728292A2B2C2D2E2F303132333435363738393A3B3C3D3E3F40414243444546474849"
"4A4B4C4D4E4F505152535455565758595A5B5C5D5E5F606162636465666768696A6B6C6D6E6F707172737475767778797A7B7C7D7E7F808182838485868788898A8B8C8D8E8F90919293"
"9495969798999A9B9C9D9E9FA0A1A2A3A4A5A6A7A8A9AAABACADAEAFB0B1B2B3B4B5B6B7B8B9BABBBCBDBEBFC0C1C2C3C4C5C6C7C8C9CACBCCCDCECFD0D1D2D3D4D5D6D7D8D9DADBDCDD"
"DEDFE0E1E2E3E4E5E6E7E8E9EAEBECEDEEEFF0F1F2F3F4F5F6F7F8F9FAFBFCFDFEFF";

#if defined(SKIP_48_ROWS) || defined(SKIP_26_ROWS) || defined(SKIP_7_ROWS) || \
    defined(SKIP_48) || defined(SKIP_26) || defined(SKIP_7)
#error SKIP_* already defined
#else
#define SKIP_7_ROWS {}, {}, {}, {}, {}, {}, {}
#define SKIP_26_ROWS SKIP_7_ROWS, SKIP_7_ROWS, SKIP_7_ROWS, {}, {}, {}, {}, {}
#define SKIP_48_ROWS SKIP_26_ROWS, SKIP_7_ROWS, SKIP_7_ROWS, SKIP_7_ROWS, {}
#define SKIP_7 0, 0, 0, 0, 0, 0, 0
#define SKIP_26 SKIP_7, SKIP_7, SKIP_7, 0, 0, 0, 0, 0
#define SKIP_48 SKIP_26, SKIP_7, SKIP_7, SKIP_7, 0
static uint8_t const kHexadecimalDigitsToByte[256][256] = {
    SKIP_48_ROWS,
    {SKIP_48, '\x00', '\x01', '\x02', '\x03', '\x04', '\x05', '\x06', '\x07',
    '\x08', '\x09', SKIP_7, '\x0a', '\x0b', '\x0c', '\x0d', '\x0e', '\x0f',
    SKIP_26, '\x0a', '\x0b', '\x0c', '\x0d', '\x0e', '\x0f'},
    {SKIP_48, '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17',
    '\x18', '\x19', SKIP_7, '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f',
    SKIP_26, '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f'},
    {SKIP_48, '\x20', '\x21', '\x22', '\x23', '\x24', '\x25', '\x26', '\x27',
    '\x28', '\x29', SKIP_7, '\x2a', '\x2b', '\x2c', '\x2d', '\x2e', '\x2f',
    SKIP_26, '\x2a', '\x2b', '\x2c', '\x2d', '\x2e', '\x2f'},
    {SKIP_48, '\x30', '\x31', '\x32', '\x33', '\x34', '\x35', '\x36', '\x37',
    '\x38', '\x39', SKIP_7, '\x3a', '\x3b', '\x3c', '\x3d', '\x3e', '\x3f',
    SKIP_26, '\x3a', '\x3b', '\x3c', '\x3d', '\x3e', '\x3f'},
    {SKIP_48, '\x40', '\x41', '\x42', '\x43', '\x44', '\x45', '\x46', '\x47',
    '\x48', '\x49', SKIP_7, '\x4a', '\x4b', '\x4c', '\x4d', '\x4e', '\x4f',
    SKIP_26, '\x4a', '\x4b', '\x4c', '\x4d', '\x4e', '\x4f'},
    {SKIP_48, '\x50', '\x51', '\x52', '\x53', '\x54', '\x55', '\x56', '\x57',
    '\x58', '\x59', SKIP_7, '\x5a', '\x5b', '\x5c', '\x5d', '\x5e', '\x5f',
    SKIP_26, '\x5a', '\x5b', '\x5c', '\x5d', '\x5e', '\x5f'},
    {SKIP_48, '\x60', '\x61', '\x62', '\x63', '\x64', '\x65', '\x66', '\x67',
    '\x68', '\x69', SKIP_7, '\x6a', '\x6b', '\x6c', '\x6d', '\x6e', '\x6f',
    SKIP_26, '\x6a', '\x6b', '\x6c', '\x6d', '\x6e', '\x6f'},
    {SKIP_48, '\x70', '\x71', '\x72', '\x73', '\x74', '\x75', '\x76', '\x77',
    '\x78', '\x79', SKIP_7, '\x7a', '\x7b', '\x7c', '\x7d', '\x7e', '\x7f',
    SKIP_26, '\x7a', '\x7b', '\x7c', '\x7d', '\x7e', '\x7f'},
    {SKIP_48, '\x80', '\x81', '\x82', '\x83', '\x84', '\x85', '\x86', '\x87',
    '\x88', '\x89', SKIP_7, '\x8a', '\x8b', '\x8c', '\x8d', '\x8e', '\x8f',
    SKIP_26, '\x8a', '\x8b', '\x8c', '\x8d', '\x8e', '\x8f'},
    {SKIP_48, '\x90', '\x91', '\x92', '\x93', '\x94', '\x95', '\x96', '\x97',
    '\x98', '\x99', SKIP_7, '\x9a', '\x9b', '\x9c', '\x9d', '\x9e', '\x9f',
    SKIP_26, '\x9a', '\x9b', '\x9c', '\x9d', '\x9e', '\x9f'},
    SKIP_7_ROWS,
    {SKIP_48, '\xa0', '\xa1', '\xa2', '\xa3', '\xa4', '\xa5', '\xa6', '\xa7',
    '\xa8', '\xa9', SKIP_7, '\xaa', '\xab', '\xac', '\xad', '\xae', '\xaf',
    SKIP_26, '\xaa', '\xab', '\xac', '\xad', '\xae', '\xaf'},
    {SKIP_48, '\xb0', '\xb1', '\xb2', '\xb3', '\xb4', '\xb5', '\xb6', '\xb7',
    '\xb8', '\xb9', SKIP_7, '\xba', '\xbb', '\xbc', '\xbd', '\xbe', '\xbf',
    SKIP_26, '\xba', '\xbb', '\xbc', '\xbd', '\xbe', '\xbf'},
    {SKIP_48, '\xc0', '\xc1', '\xc2', '\xc3', '\xc4', '\xc5', '\xc6', '\xc7',
    '\xc8', '\xc9', SKIP_7, '\xca', '\xcb', '\xcc', '\xcd', '\xce', '\xcf',
    SKIP_26, '\xca', '\xcb', '\xcc', '\xcd', '\xce', '\xcf'},
    {SKIP_48, '\xd0', '\xd1', '\xd2', '\xd3', '\xd4', '\xd5', '\xd6', '\xd7',
    '\xd8', '\xd9', SKIP_7, '\xda', '\xdb', '\xdc', '\xdd', '\xde', '\xdf',
    SKIP_26, '\xda', '\xdb', '\xdc', '\xdd', '\xde', '\xdf'},
    {SKIP_48, '\xe0', '\xe1', '\xe2', '\xe3', '\xe4', '\xe5', '\xe6', '\xe7',
    '\xe8', '\xe9', SKIP_7, '\xea', '\xeb', '\xec', '\xed', '\xee', '\xef',
    SKIP_26, '\xea', '\xeb', '\xec', '\xed', '\xee', '\xef'},
    {SKIP_48, '\xf0', '\xf1', '\xf2', '\xf3', '\xf4', '\xf5', '\xf6', '\xf7',
    '\xf8', '\xf9', SKIP_7, '\xfa', '\xfb', '\xfc', '\xfd', '\xfe', '\xff',
    SKIP_26, '\xfa', '\xfb', '\xfc', '\xfd', '\xfe', '\xff'},
    SKIP_26_ROWS,
    {SKIP_48, '\xa0', '\xa1', '\xa2', '\xa3', '\xa4', '\xa5', '\xa6', '\xa7',
    '\xa8', '\xa9', SKIP_7, '\xaa', '\xab', '\xac', '\xad', '\xae', '\xaf',
    SKIP_26, '\xaa', '\xab', '\xac', '\xad', '\xae', '\xaf'},
    {SKIP_48, '\xb0', '\xb1', '\xb2', '\xb3', '\xb4', '\xb5', '\xb6', '\xb7',
    '\xb8', '\xb9', SKIP_7, '\xba', '\xbb', '\xbc', '\xbd', '\xbe', '\xbf',
    SKIP_26, '\xba', '\xbb', '\xbc', '\xbd', '\xbe', '\xbf'},
    {SKIP_48, '\xc0', '\xc1', '\xc2', '\xc3', '\xc4', '\xc5', '\xc6', '\xc7',
    '\xc8', '\xc9', SKIP_7, '\xca', '\xcb', '\xcc', '\xcd', '\xce', '\xcf',
    SKIP_26, '\xca', '\xcb', '\xcc', '\xcd', '\xce', '\xcf'},
    {SKIP_48, '\xd0', '\xd1', '\xd2', '\xd3', '\xd4', '\xd5', '\xd6', '\xd7',
    '\xd8', '\xd9', SKIP_7, '\xda', '\xdb', '\xdc', '\xdd', '\xde', '\xdf',
    SKIP_26, '\xda', '\xdb', '\xdc', '\xdd', '\xde', '\xdf'},
    {SKIP_48, '\xe0', '\xe1', '\xe2', '\xe3', '\xe4', '\xe5', '\xe6', '\xe7',
    '\xe8', '\xe9', SKIP_7, '\xea', '\xeb', '\xec', '\xed', '\xee', '\xef',
    SKIP_26, '\xea', '\xeb', '\xec', '\xed', '\xee', '\xef'},
    {SKIP_48, '\xf0', '\xf1', '\xf2', '\xf3', '\xf4', '\xf5', '\xf6', '\xf7',
    '\xf8', '\xf9', SKIP_7, '\xfa', '\xfb', '\xfc', '\xfd', '\xfe', '\xff',
    SKIP_26, '\xfa', '\xfb', '\xfc', '\xfd', '\xfe', '\xff'}};

static uint8_t const kHexadecimalDigitsToNibble[256] = {
    SKIP_48, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    SKIP_7, '\xa', '\xb', '\xc', '\xd', '\xe', '\xf',
    SKIP_26, '\xa', '\xb', '\xc', '\xd', '\xe', '\xf'};
#undef SKIP_7_ROWS
#undef SKIP_26_ROWS
#undef SKIP_48_ROWS
#undef SKIP_7
#undef SKIP_26
#undef SKIP_48
#endif


template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value,
    void>
HexadecimalEncode(Container const& input, not_null<Container*> output) {
  Container const& bytes = input;
  Container digits;
  digits.resize(bytes.size() << 1);
  // The following was undefined behaviour pre-C++11, but it is now well-defined
  // even when |digits.size() == 0|.  We do not use |digits.data()| because this
  // only works for |std::vector| (it is read-only in a |std::basic_string|).
  char* digit = &digits[0];
  for (uint8_t const byte : bytes) {
    // The following is four times faster than copying both bytes by hand.
    memcpy(digit, &kByteToHexadecimalDigits[byte << 1], 2);
    digit += 2;
  }
  *output = std::move(digits);
}

template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value>
HexadecimalDecode(Container const& input, not_null<Container*> output) {;
  Container const& digits = input;
  Container& bytes = *output;
  std::size_t const odd = digits.size() % 2;
  CHECK(bytes.data() >= &digits.data()[input.size()] ||
        bytes.data() <= &digits.data()[1]);
  // Do not shrink |bytes| in case we are decoding in-place.
  bytes.resize(std::max(bytes.size(), digits.size() / 2));
  auto byte = bytes.begin();
  for (auto digit = digits.begin();
       digit != digits.end() - odd;
       ++digit, ++byte) {
    *byte = (kHexadecimalDigitsToNibble[*digit] << 4) + 
             kHexadecimalDigitsToNibble[*++digit];
  }
  bytes.resize(digits.size() / 2);
}

}  // namespace base
}  // namespace principia
