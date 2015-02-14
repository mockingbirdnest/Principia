
#include "base/hexadecimal.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

class HexadecimalTest : public testing::Test {
 protected:
  std::string const bytes_ = std::string("\0\x7F\x80\xFFgh\n\7", 8);
  std::string const lowercase_digits_ = "00""7f""80""ff""67""68""0a""07";
  std::string const uppercase_digits_ = "00""7F""80""FF""67""68""0A""07";
};

TEST_F(HexadecimalTest, EncodeAndDecode) {
  std::string digits;
  HexadecimalEncode<std::string>(bytes_, &digits);
  EXPECT_EQ(uppercase_digits_, digits);
  std::string bytes;
  HexadecimalDecode<std::string>(digits, &bytes);
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  std::string str = bytes_;
  HexadecimalEncode<std::string>(str, &str);
  EXPECT_EQ(uppercase_digits_, str);
  HexadecimalDecode<std::string>(str, &str);
  EXPECT_EQ(bytes_, str);
}

TEST_F(HexadecimalTest, CaseInsensitive) {
  std::string bytes;
  HexadecimalDecode<std::string>(lowercase_digits_, &bytes);
  EXPECT_EQ(bytes_, bytes);
  HexadecimalDecode<std::string>(uppercase_digits_, &bytes);
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, Invalid) {
  std::string bytes;
  std::string digits = "abc";
  HexadecimalDecode<std::string>(digits, &bytes);;
  EXPECT_EQ("\xAB", bytes);
  digits = "0agcde";
  HexadecimalDecode<std::string>(digits, &bytes);
  EXPECT_EQ(std::string("\x0A\x0C\xDE", 3), bytes);
}

}  // namespace base
}  // namespace principia
