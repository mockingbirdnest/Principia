
#include "base/hexadecimal.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::ElementsAre;

namespace principia {
namespace base {

class HexadecimalTest : public testing::Test {
 protected:
  HexadecimalTest() {
    std::string str = std::string("\0\x7F\x80\xFFgh\n\7", kBytes);
    bytes_.assign(str.begin(), str.end());
    str = "00""7f""80""ff""67""68""0a""07";
    lowercase_digits_.assign(str.begin(), str.end());
    str = "00""7F""80""FF""67""68""0A""07";
    uppercase_digits_.assign(str.begin(), str.end());
  }
  static std::size_t const kBytes = 8;
  static std::size_t const kDigits = kBytes << 1;
  std::vector<uint8_t> bytes_;
  std::vector<uint8_t> lowercase_digits_;
  std::vector<uint8_t> uppercase_digits_;
};

using HexadecimalDeathTest = HexadecimalTest;

TEST_F(HexadecimalTest, EncodeAndDecode) {
  std::vector<uint8_t> digits(kDigits);
  HexadecimalEncode(&bytes_[0], bytes_.size(), &digits[0], digits.size());
  EXPECT_EQ(uppercase_digits_, digits);
  std::vector<uint8_t> bytes(kBytes);
  HexadecimalDecode(&digits[0], digits.size(), &bytes[0], bytes.size());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  std::vector<uint8_t> str(kDigits);
  memcpy(&str[kBytes], &bytes_[0], kBytes);
  HexadecimalEncode(&str[kBytes], kBytes, &str[0], str.size());
  EXPECT_EQ(uppercase_digits_, str);
  HexadecimalDecode(&str[0], str.size(), &str[1], kBytes);
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&str[1], &str[kBytes + 1]));
  str = uppercase_digits_;
  HexadecimalDecode(&str[0], str.size(), &str[0], kBytes);
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&str[0], &str[kBytes]));
}

TEST_F(HexadecimalTest, CaseInsensitive) {
  std::vector<uint8_t> bytes(kBytes);
  HexadecimalDecode(&lowercase_digits_[0], lowercase_digits_.size(),
                    &bytes[0], bytes.size());
  EXPECT_EQ(bytes_, bytes);
  HexadecimalDecode(&uppercase_digits_[0], uppercase_digits_.size(),
                    &bytes[0], bytes.size());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, Invalid) {
  std::vector<uint8_t> bytes(1);
  std::vector<uint8_t> digits = {'a', 'b', 'c'};
  HexadecimalDecode(&digits[0], digits.size(), &bytes[0], bytes.size());
  EXPECT_THAT(bytes, ElementsAre('\xAB'));;
  digits = {'0', 'a', 'g', 'c', 'd', 'e'};
  bytes.resize(3);
  HexadecimalDecode(&digits[0], digits.size(), &bytes[0], bytes.size());
  EXPECT_THAT(bytes, ElementsAre('\x0A', '\x0C', '\xDE'));
}

}  // namespace base
}  // namespace principia
