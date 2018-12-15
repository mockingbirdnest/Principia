
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
    : bytes_("\x00\x7F\x80\xFF\x67\x68\x0A\x07"),
      lowercase_digits_("00""7f""80""ff""67""68""0a""07"),
      uppercase_digits_("00""7F""80""FF""67""68""0A""07"),
      digits_(digit_count) {
    CHECK_EQ(bytes_.size, byte_count);
    CHECK_EQ(lowercase_digits_.size, digit_count);
    CHECK_EQ(uppercase_digits_.size, digit_count);
  }

  static std::int64_t const byte_count = 8;
  static std::int64_t const digit_count = byte_count << 1;

  Array<std::uint8_t const> const bytes_;
  Array<char const> const lowercase_digits_;
  Array<char const> const uppercase_digits_;
  UniqueArray<char> const digits_;
  HexadecimalEncoder</*null_terminated=*/false> encoder_;
};

using HexadecimalDeathTest = HexadecimalTest;

TEST_F(HexadecimalTest, EncodeAndDecode) {
  encoder_.Encode(bytes_, digits_.get());
  EXPECT_EQ(uppercase_digits_, digits_);
  EXPECT_EQ(uppercase_digits_.size, encoder_.EncodedLength(bytes_));
  UniqueArray<std::uint8_t> bytes(byte_count);
  encoder_.Decode(digits_.get(), bytes.get());
  EXPECT_EQ(bytes_, bytes.get());
  EXPECT_EQ(bytes_.size, encoder_.DecodedLength(digits_.get()));
}

TEST_F(HexadecimalTest, UniqueEncodeAndDecode) {
  auto const digits = encoder_.Encode(bytes_);
  EXPECT_EQ(uppercase_digits_, digits);
  auto const bytes = encoder_.Decode(digits.get());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  auto const buffer = std::make_unique<std::uint8_t[]>(digit_count);
  auto const buffer_characters = reinterpret_cast<char*>(buffer.get());
  std::memcpy(&buffer[1], bytes_.data, byte_count);
  encoder_.Encode({&buffer[1], byte_count},
                  {&buffer_characters[0], digit_count});
  EXPECT_EQ(uppercase_digits_, Array<char>(&buffer_characters[0], digit_count));
  std::memcpy(&buffer[0], bytes_.data, byte_count);
  encoder_.Encode({&buffer[0], byte_count},
                  {&buffer_characters[0], digit_count});
  EXPECT_EQ(uppercase_digits_, Array<std::uint8_t>(&buffer[0], digit_count));
  encoder_.Decode({&buffer_characters[0], digit_count},
                  {&buffer[1], byte_count});
  EXPECT_EQ(bytes_, Array<std::uint8_t>(&buffer[1], byte_count));
  std::memcpy(&buffer[0], uppercase_digits_.data, digit_count);
  encoder_.Decode({&buffer_characters[0], digit_count},
                  {&buffer[0], byte_count});
  EXPECT_EQ(bytes_, Array<std::uint8_t>(&buffer[0], byte_count));
}

TEST_F(HexadecimalTest, LargeOutput) {
  std::int64_t const digits_size = digit_count + 42;
  UniqueArray<char> const digits(digits_size);
  std::memset(digits.data.get(), 'X', digits_size);
  encoder_.Encode(bytes_, digits.get());
  EXPECT_EQ(uppercase_digits_, Array<char>(digits.data.get(), digit_count));
  EXPECT_THAT(std::vector<std::uint8_t>(&digits.data[digit_count],
                                        &digits.data[digits_size]),
              Each('X'));
  std::int64_t const bytes_size = byte_count + 42;
  UniqueArray<std::uint8_t> bytes(bytes_size);
  std::memset(bytes.data.get(), 'Y', bytes_size);
  encoder_.Decode(uppercase_digits_, bytes.get());
  EXPECT_EQ(bytes_, Array<std::uint8_t>(bytes.data.get(), byte_count));
  EXPECT_THAT(std::vector<std::uint8_t>(&bytes.data[byte_count],
                                        &bytes.data[bytes_size]),
              Each('Y'));
}

TEST_F(HexadecimalTest, Adjacent) {
  auto const buffer =
      std::make_unique<std::uint8_t[]>(digit_count + byte_count);
  auto const buffer_characters = reinterpret_cast<char*>(buffer.get());
  std::memcpy(&buffer[0], bytes_.data, byte_count);
  encoder_.Encode({&buffer[0], byte_count},
                  {&buffer_characters[byte_count], digit_count});
  EXPECT_EQ(uppercase_digits_,
            Array<std::uint8_t>(&buffer[byte_count], digit_count));
  std::memcpy(&buffer[0], uppercase_digits_.data, digit_count);
  encoder_.Decode({&buffer_characters[0], digit_count},
                  {&buffer[digit_count], byte_count});
  EXPECT_EQ(bytes_, Array<std::uint8_t>(&buffer[digit_count], byte_count));
  encoder_.Decode({&buffer_characters[0], digit_count + 1},
                  {&buffer[digit_count], byte_count});
  EXPECT_EQ(bytes_, Array<std::uint8_t>(&buffer[digit_count], byte_count));
}

TEST_F(HexadecimalDeathTest, Overlap) {
  auto const buffer =
      std::make_unique<std::uint8_t[]>(digit_count + byte_count - 1);
  auto const buffer_characters = reinterpret_cast<char*>(buffer.get());
  EXPECT_DEATH({
    encoder_.Encode({&buffer[digit_count - 1], byte_count},
                    {&buffer_characters[0], digit_count});
  }, "bad overlap");
  EXPECT_DEATH({
    encoder_.Decode({&buffer_characters[0], digit_count},
                    {&buffer[digit_count - 1], byte_count});
  }, "bad overlap");
  EXPECT_DEATH({
    encoder_.Encode({&buffer[2], byte_count},
                    {&buffer_characters[0], digit_count});
  }, "bad overlap");
  EXPECT_DEATH({
    encoder_.Decode({&buffer_characters[0], digit_count},
                    {&buffer[2], byte_count});
  }, "bad overlap");
}

TEST_F(HexadecimalDeathTest, Size) {
  std::vector<std::uint8_t> bytes(byte_count);
  std::vector<char> digits(digit_count);
  EXPECT_DEATH({
    encoder_.Encode(bytes, {digits.data(), digits.size() - 1});
  }, "too small");
  EXPECT_DEATH({
    encoder_.Decode(digits, {bytes.data(), bytes.size() - 1});
  }, "too small");
}

TEST_F(HexadecimalTest, CaseInsensitive) {
  std::vector<std::uint8_t> bytes(byte_count);
  encoder_.Decode(lowercase_digits_, bytes);
  EXPECT_EQ(bytes_, Array<std::uint8_t>(bytes));
  encoder_.Decode(uppercase_digits_, bytes);
  EXPECT_EQ(bytes_, Array<std::uint8_t>(bytes));
}

TEST_F(HexadecimalTest, Invalid) {
  std::vector<std::uint8_t> bytes(1);
  std::vector<char> digits = {'a', 'b', 'c'};
  encoder_.Decode(digits, bytes);
  EXPECT_THAT(bytes, ElementsAre('\xAB'));
  digits = {'0', 'a', 'g', 'c', 'd', 'e'};
  bytes.resize(3);
  encoder_.Decode(digits, bytes);
  EXPECT_THAT(bytes, ElementsAre('\x0A', '\x0C', '\xDE'));
}

}  // namespace base
}  // namespace principia
