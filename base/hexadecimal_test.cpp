
#include "base/hexadecimal.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Each;
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

  static int64_t const kBytes = 8;
  static int64_t const kDigits = kBytes << 1;

  std::vector<uint8_t> bytes_;
  std::vector<uint8_t> lowercase_digits_;
  std::vector<uint8_t> uppercase_digits_;
};

using HexadecimalDeathTest = HexadecimalTest;

TEST_F(HexadecimalTest, EncodeAndDecode) {
  std::vector<uint8_t> digits(kDigits);
  HexadecimalEncode(bytes_.data(), bytes_.size(), digits.data(), digits.size());
  EXPECT_EQ(uppercase_digits_, digits);
  std::vector<uint8_t> bytes(kBytes);
  HexadecimalDecode(digits.data(), digits.size(), bytes.data(), bytes.size());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  std::vector<uint8_t> buffer(kDigits);
  memcpy(&buffer[1], bytes_.data(), kBytes);
  HexadecimalEncode(&buffer[1], kBytes, &buffer[0], kDigits);
  EXPECT_EQ(uppercase_digits_, buffer);
  memcpy(&buffer[0], bytes_.data(), kBytes);
  HexadecimalEncode(&buffer[0], kBytes, &buffer[0], kDigits);
  EXPECT_EQ(uppercase_digits_, buffer);
  HexadecimalDecode(&buffer[0], kDigits, &buffer[1], kBytes);
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&buffer[1], &buffer[kBytes + 1]));
  buffer = uppercase_digits_;
  HexadecimalDecode(&buffer[0], kDigits, &buffer[0], kBytes);
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&buffer[0], &buffer[kBytes]));
}

TEST_F(HexadecimalTest, LargeOutput) {
  std::vector<uint8_t> digits(kDigits + 42, 'X');
  HexadecimalEncode(bytes_.data(), bytes_.size(), digits.data(), digits.size());
  EXPECT_EQ(uppercase_digits_,
            std::vector<uint8_t>(&digits[0], &digits[kDigits]));
  EXPECT_THAT(std::vector<uint8_t>(&digits[kDigits], &digits[digits.size()]),
              Each('X'));
  std::vector<uint8_t> bytes(kBytes + 42, 'Y');
  HexadecimalDecode(uppercase_digits_.data(), uppercase_digits_.size(),
                    bytes.data(), bytes.size());
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&bytes[0], &bytes[kBytes]));
  EXPECT_THAT(std::vector<uint8_t>(&bytes[kBytes], &bytes[bytes.size()]),
              Each('Y'));
}

TEST_F(HexadecimalTest, Adjacent) {
  std::vector<uint8_t> buffer(kDigits + kBytes);
  memcpy(&buffer[0], bytes_.data(), kBytes);
  HexadecimalEncode(&buffer[0], kBytes, &buffer[kBytes], kDigits);
  EXPECT_EQ(uppercase_digits_,
            std::vector<uint8_t>(&buffer[kBytes], &buffer[kBytes + kDigits]));
  memcpy(&buffer[0], uppercase_digits_.data(), kDigits);
  HexadecimalDecode(&buffer[0], kDigits, &buffer[kDigits], kBytes);
  EXPECT_EQ(bytes_,
            std::vector<uint8_t>(&buffer[kDigits], &buffer[kDigits + kBytes]));
  HexadecimalDecode(&buffer[0], kDigits + 1, &buffer[kDigits], kBytes);
  EXPECT_EQ(bytes_,
            std::vector<uint8_t>(&buffer[kDigits], &buffer[kDigits + kBytes]));
}

TEST_F(HexadecimalDeathTest, Overlap) {
  std::vector<uint8_t> buffer(kDigits + kBytes - 1);
  EXPECT_DEATH({
    HexadecimalEncode(&buffer[kDigits - 1], kBytes, &buffer[0], kDigits);
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode(&buffer[0], kDigits, &buffer[kDigits - 1], kBytes);
  }, "bad overlap");
  buffer = std::vector<uint8_t>(kDigits);
  EXPECT_DEATH({
    HexadecimalEncode(&buffer[2], kBytes, &buffer[0], kDigits);
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode(&buffer[0], kDigits, &buffer[2], kBytes);
  }, "bad overlap");
}

TEST_F(HexadecimalDeathTest, Size) {
  std::vector<uint8_t> bytes(kBytes);
  std::vector<uint8_t> digits(kDigits);
  EXPECT_DEATH({
    HexadecimalEncode(bytes.data(), bytes.size(),
                      digits.data(), digits.size() - 1);
  }, "too small");
  EXPECT_DEATH({
    HexadecimalDecode(digits.data(), digits.size(),
                      bytes.data(), bytes.size() - 1);
  }, "too small");
}

TEST_F(HexadecimalTest, CaseInsensitive) {
  std::vector<uint8_t> bytes(kBytes);
  HexadecimalDecode(lowercase_digits_.data(), lowercase_digits_.size(),
                    bytes.data(), bytes.size());
  EXPECT_EQ(bytes_, bytes);
  HexadecimalDecode(uppercase_digits_.data(), uppercase_digits_.size(),
                    bytes.data(), bytes.size());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, Invalid) {
  std::vector<uint8_t> bytes(1);
  std::vector<uint8_t> digits = {'a', 'b', 'c'};
  HexadecimalDecode(digits.data(), digits.size(), bytes.data(), bytes.size());
  EXPECT_THAT(bytes, ElementsAre('\xAB'));
  digits = {'0', 'a', 'g', 'c', 'd', 'e'};
  bytes.resize(3);
  HexadecimalDecode(digits.data(), digits.size(), bytes.data(), bytes.size());
  EXPECT_THAT(bytes, ElementsAre('\x0A', '\x0C', '\xDE'));
}

}  // namespace base
}  // namespace principia
