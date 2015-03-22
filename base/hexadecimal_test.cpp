
#include "base/hexadecimal.hpp"

#include <memory>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Each;
using testing::ElementsAre;

namespace principia {
namespace base {

class HexadecimalTest : public testing::Test {
 protected:
  HexadecimalTest() 
    : bytes_(Bytes::New(std::string("\0\x7F\x80\xFFgh\n\7", kBytes))),
      digits_(Bytes::New(kDigits)),
      lowercase_digits_(Bytes::New("00""7f""80""ff""67""68""0a""07")),
      uppercase_digits_(Bytes::New("00""7F""80""FF""67""68""0A""07")) {}

  static int64_t const kBytes = 8;
  static int64_t const kDigits = kBytes << 1;

  UniqueBytes bytes_;
  UniqueBytes digits_;
  UniqueBytes lowercase_digits_;
  UniqueBytes uppercase_digits_;
};

using HexadecimalDeathTest = HexadecimalTest;

TEST_F(HexadecimalTest, EncodeAndDecode) {
  HexadecimalEncode(*bytes_, digits_.get());
  EXPECT_EQ(uppercase_digits_, digits_);
  UniqueBytes bytes(Bytes::New(kBytes));
  HexadecimalDecode(*digits_, bytes.get());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  auto buffer = std::make_unique<uint8_t[]>(kDigits);
  memcpy(&buffer[1], bytes_.data, kBytes);
  HexadecimalEncode({&buffer[1], kBytes}, {&buffer[0], kDigits});
  EXPECT_EQ(uppercase_digits_,
            std::vector<uint8_t>(&buffer[0], &buffer[kDigits]));
  memcpy(&buffer[0], bytes_.data, kBytes);
  HexadecimalEncode({&buffer[0], kBytes}, {&buffer[0], kDigits});
  EXPECT_EQ(uppercase_digits_,
            std::vector<uint8_t>(&buffer[0], &buffer[kDigits]));
  HexadecimalDecode({&buffer[0], kDigits}, {&buffer[1], kBytes});
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&buffer[1], &buffer[kBytes + 1]));
  memcpy(&buffer[0], uppercase_digits_.data, kDigits);
  HexadecimalDecode({&buffer[0], kDigits}, {&buffer[0], kBytes});
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&buffer[0], &buffer[kBytes]));
}

TEST_F(HexadecimalTest, LargeOutput) {
  int64_t const digits_size = kDigits + 42;
  auto digits = std::make_unique<uint8_t[]>(digits_size);
  memset(&digits[0], 'X', digits_size);
  HexadecimalEncode(bytes_, &digits);
  EXPECT_EQ(uppercase_digits_,
            std::vector<uint8_t>(&digits[0], &digits[kDigits]));
  EXPECT_THAT(std::vector<uint8_t>(&digits[kDigits], &digits[digits_size]),
              Each('X'));
  int64_t const bytes_size = kBytes + 42;
  auto bytes = std::make_unique<uint8_t[]>(bytes_size);
  memset(&bytes[0], 'Y', bytes_size);
  HexadecimalDecode(uppercase_digits_.data(), uppercase_digits_.size(),
                    &bytes[0], bytes_size);
  EXPECT_EQ(bytes_, std::vector<uint8_t>(&bytes[0], &bytes[kBytes]));
  EXPECT_THAT(std::vector<uint8_t>(&bytes[kBytes], &bytes[bytes_size]),
              Each('Y'));
}

TEST_F(HexadecimalTest, Adjacent) {
  auto buffer = std::make_unique<uint8_t[]>(kDigits + kBytes);
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
  auto buffer = std::make_unique<uint8_t[]>(kDigits + kBytes - 1);
  EXPECT_DEATH({
    HexadecimalEncode(&buffer[kDigits - 1], kBytes, &buffer[0], kDigits);
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode(&buffer[0], kDigits, &buffer[kDigits - 1], kBytes);
  }, "bad overlap");
  buffer = std::make_unique<uint8_t[]>(kDigits);
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
