#pragma once

#include <string.h>

#include "base/hexadecimal.hpp"

namespace principia {
namespace base {

static char const* kByteToHexadecimalDigits[256] = {
    "00", "01", "02", "03", "04", "05", "06", "07",
    "08", "09", "0A", "0B", "0C", "0D", "0E", "0F",
    "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "1A", "1B", "1C", "1D", "1E", "1F",
    "20", "21", "22", "23", "24", "25", "26", "27",
    "28", "29", "2A", "2B", "2C", "2D", "2E", "2F",
    "30", "31", "32", "33", "34", "35", "36", "37",
    "38", "39", "3A", "3B", "3C", "3D", "3E", "3F",
    "40", "41", "42", "43", "44", "45", "46", "47",
    "48", "49", "4A", "4B", "4C", "4D", "4E", "4F",
    "50", "51", "52", "53", "54", "55", "56", "57",
    "58", "59", "5A", "5B", "5C", "5D", "5E", "5F",
    "60", "61", "62", "63", "64", "65", "66", "67",
    "68", "69", "6A", "6B", "6C", "6D", "6E", "6F",
    "70", "71", "72", "73", "74", "75", "76", "77",
    "78", "79", "7A", "7B", "7C", "7D", "7E", "7F",
    "80", "81", "82", "83", "84", "85", "86", "87",
    "88", "89", "8A", "8B", "8C", "8D", "8E", "8F",
    "90", "91", "92", "93", "94", "95", "96", "97",
    "98", "99", "9A", "9B", "9C", "9D", "9E", "9F",
    "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7",
    "A8", "A9", "AA", "AB", "AC", "AD", "AE", "AF",
    "B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7",
    "B8", "B9", "BA", "BB", "BC", "BD", "BE", "BF",
    "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7",
    "C8", "C9", "CA", "CB", "CC", "CD", "CE", "CF",
    "D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7",
    "D8", "D9", "DA", "DB", "DC", "DD", "DE", "DF",
    "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7",
    "E8", "E9", "EA", "EB", "EC", "ED", "EE", "EF",
    "F0", "F1", "F2", "F3", "F4", "F5", "F6", "F7",
    "F8", "F9", "FA", "FB", "FC", "FD", "FE", "FF"};

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
  digits.resize(bytes.size() * 2);
  // The following was undefined behaviour pre-C++11, but it is now well-defined
  // even when |digits.size() == 0|.  We do not use |digits.data()| because this
  // only works for |std::vector| (it is read-only in a |std::basic_string|).
  char* digit = &digits[0];
  for (uint8_t const byte : bytes) {
    // The following is four times faster than copying both bytes by hand.
    memcpy(digit, kByteToHexadecimalDigits[byte], 2);
    digit += 2;
  }
  *output = std::move(digits);
}

template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value,
    bool>
HexadecimalDecode(Container const& input, not_null<Container*> output) {
  bool valid = true;
  Container const& digits = input;
  Container& bytes = *output;
  if (digits.size() % 2 != 0) {
    valid = false;
  }
  // Do not shrink |bytes| in case we are decoding in-place.
  bytes.resize(std::max(bytes.size(), digits.size() / 2));
  auto byte = bytes.begin();
  for (auto digit = digits.begin();
       digit != digits.end() - digits.size() % 2;
       ++digit, ++byte) {
    *byte = kHexadecimalDigitsToByte[*digit][*++digit];
    if (*byte == 0 && (*digit != '0' || *(digit - 1) != '0')) {
      valid = false;
    }
  }
  bytes.resize(digits.size() / 2);
  return valid;
}

}  // namespace base
}  // namespace principia
