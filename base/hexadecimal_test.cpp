
#include "base/hexadecimal.hpp"

#include <memory>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Each;
using testing::ElementsAre;

namespace principia {
namespace base {

class HexadecimalTest : public testing::Test {
 protected:
  HexadecimalTest()
    : bytes_(byte_count),
      digits_(digit_count) {
    std::string const lowercase_digits = "00""7f""80""ff""67""68""0a""07";
    std::string const uppercase_digits = "00""7F""80""FF""67""68""0A""07";
    std::memcpy(bytes_.data.get(), "\0\x7F\x80\xFFgh\n\7", byte_count);
    lowercase_digits_ = UniqueBytes(lowercase_digits.size());
    uppercase_digits_ = UniqueBytes(uppercase_digits.size());
    std::memcpy(lowercase_digits_.data.get(),
                lowercase_digits.c_str(),
                lowercase_digits.size());
    std::memcpy(uppercase_digits_.data.get(),
                uppercase_digits.c_str(),
                uppercase_digits.size());
  }

  static std::int64_t const byte_count = 8;
  static std::int64_t const digit_count = byte_count << 1;

  UniqueBytes bytes_;
  UniqueBytes digits_;
  UniqueBytes lowercase_digits_;
  UniqueBytes uppercase_digits_;
};

using HexadecimalDeathTest = HexadecimalTest;

TEST_F(HexadecimalTest, EncodeAndDecode) {
  HexadecimalEncode(bytes_.get(), digits_.get());
  EXPECT_EQ(uppercase_digits_, digits_);
  UniqueBytes bytes(byte_count);
  HexadecimalDecode(digits_.get(), bytes.get());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, UniqueEncodeAndDecode) {
  auto const digits =
      HexadecimalEncode(bytes_.get(), /*null_terminated=*/false);
  EXPECT_EQ(uppercase_digits_, digits);
  auto const bytes = HexadecimalDecode(digits.get());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  auto buffer = std::make_unique<std::uint8_t[]>(digit_count);
  std::memcpy(&buffer[1], bytes_.data.get(), byte_count);
  HexadecimalEncode({&buffer[1], byte_count}, {&buffer[0], digit_count});
  EXPECT_EQ(uppercase_digits_, Bytes(&buffer[0], digit_count));
  std::memcpy(&buffer[0], bytes_.data.get(), byte_count);
  HexadecimalEncode({&buffer[0], byte_count}, {&buffer[0], digit_count});
  EXPECT_EQ(uppercase_digits_, Bytes(&buffer[0], digit_count));
  HexadecimalDecode({&buffer[0], digit_count}, {&buffer[1], byte_count});
  EXPECT_EQ(bytes_, Bytes(&buffer[1], byte_count));
  std::memcpy(&buffer[0], uppercase_digits_.data.get(), digit_count);
  HexadecimalDecode({&buffer[0], digit_count}, {&buffer[0], byte_count});
  EXPECT_EQ(bytes_, Bytes(&buffer[0], byte_count));
}

TEST_F(HexadecimalTest, LargeOutput) {
  std::int64_t const digits_size = digit_count + 42;
  UniqueBytes digits(digits_size);
  std::memset(digits.data.get(), 'X', digits_size);
  HexadecimalEncode(bytes_.get(), digits.get());
  EXPECT_EQ(uppercase_digits_, Bytes(digits.data.get(), digit_count));
  EXPECT_THAT(std::vector<std::uint8_t>(&digits.data[digit_count],
                                        &digits.data[digits_size]),
              Each('X'));
  std::int64_t const bytes_size = byte_count + 42;
  UniqueBytes bytes(bytes_size);
  std::memset(bytes.data.get(), 'Y', bytes_size);
  HexadecimalDecode(uppercase_digits_.get(), bytes.get());
  EXPECT_EQ(bytes_, Bytes(bytes.data.get(), byte_count));
  EXPECT_THAT(std::vector<std::uint8_t>(&bytes.data[byte_count],
                                        &bytes.data[bytes_size]),
              Each('Y'));
}

TEST_F(HexadecimalTest, Adjacent) {
  auto buffer = std::make_unique<std::uint8_t[]>(digit_count + byte_count);
  std::memcpy(&buffer[0], bytes_.data.get(), byte_count);
  HexadecimalEncode({&buffer[0], byte_count},
                    {&buffer[byte_count], digit_count});
  EXPECT_EQ(uppercase_digits_, Bytes(&buffer[byte_count], digit_count));
  std::memcpy(&buffer[0], uppercase_digits_.data.get(), digit_count);
  HexadecimalDecode({&buffer[0], digit_count},
                    {&buffer[digit_count], byte_count});
  EXPECT_EQ(bytes_, Bytes(&buffer[digit_count], byte_count));
  HexadecimalDecode({&buffer[0], digit_count + 1},
                    {&buffer[digit_count], byte_count});
  EXPECT_EQ(bytes_, Bytes(&buffer[digit_count], byte_count));
}

TEST_F(HexadecimalDeathTest, Overlap) {
  auto buffer = std::make_unique<std::uint8_t[]>(digit_count + byte_count - 1);
  EXPECT_DEATH({
    HexadecimalEncode({&buffer[digit_count - 1], byte_count},
                      {&buffer[0], digit_count});
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode({&buffer[0], digit_count},
                      {&buffer[digit_count - 1], byte_count});
  }, "bad overlap");
  buffer = std::make_unique<std::uint8_t[]>(digit_count);
  EXPECT_DEATH({
    HexadecimalEncode({&buffer[2], byte_count}, {&buffer[0], digit_count});
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode({&buffer[0], digit_count}, {&buffer[2], byte_count});
  }, "bad overlap");
}

TEST_F(HexadecimalDeathTest, Size) {
  std::vector<std::uint8_t> bytes(byte_count);
  std::vector<std::uint8_t> digits(digit_count);
  EXPECT_DEATH({
    HexadecimalEncode({bytes.data(), bytes.size()},
                      {digits.data(), digits.size() - 1});
  }, "too small");
  EXPECT_DEATH({
    HexadecimalDecode({digits.data(), digits.size()},
                      {bytes.data(), bytes.size() - 1});
  }, "too small");
}

TEST_F(HexadecimalTest, CaseInsensitive) {
  std::vector<std::uint8_t> bytes(byte_count);
  HexadecimalDecode(lowercase_digits_.get(), {bytes.data(), bytes.size()});
  EXPECT_EQ(bytes_, Bytes(bytes.data(), bytes.size()));
  HexadecimalDecode(uppercase_digits_.get(), {bytes.data(), bytes.size()});
  EXPECT_EQ(bytes_, Bytes(bytes.data(), bytes.size()));
}

TEST_F(HexadecimalTest, Invalid) {
  std::vector<std::uint8_t> bytes(1);
  std::vector<std::uint8_t> digits = {'a', 'b', 'c'};
  HexadecimalDecode({digits.data(), digits.size()},
                    {bytes.data(), bytes.size()});
  EXPECT_THAT(bytes, ElementsAre('\xAB'));
  digits = {'0', 'a', 'g', 'c', 'd', 'e'};
  bytes.resize(3);
  HexadecimalDecode({digits.data(), digits.size()},
                    {bytes.data(), bytes.size()});
  EXPECT_THAT(bytes, ElementsAre('\x0A', '\x0C', '\xDE'));
}

}  // namespace base
}  // namespace principia
