
#include "base/hexadecimal.hpp"

#include <memory>
#include <string>
#include <vector>

#include "base/bytes.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Each;
using testing::ElementsAre;

namespace principia {
namespace base {

class HexadecimalTest : public testing::Test {
 protected:
  HexadecimalTest() 
    : bytes_(kBytes),
      digits_(kDigits) {
    std::string const lowercase_digits = "00""7f""80""ff""67""68""0a""07";
    std::string const uppercase_digits = "00""7F""80""FF""67""68""0A""07";
    std::memcpy(bytes_.data.get(), "\0\x7F\x80\xFFgh\n\7", kBytes);
    lowercase_digits_ = std::move(UniqueBytes(lowercase_digits.size()));
    uppercase_digits_ = std::move(UniqueBytes(uppercase_digits.size()));
    std::memcpy(lowercase_digits_.data.get(),
                lowercase_digits.c_str(),
                lowercase_digits.size());
    std::memcpy(uppercase_digits_.data.get(),
                uppercase_digits.c_str(),
                uppercase_digits.size());
  }

  static int64_t const kBytes = 8;
  static int64_t const kDigits = kBytes << 1;

  UniqueBytes bytes_;
  UniqueBytes digits_;
  UniqueBytes lowercase_digits_;
  UniqueBytes uppercase_digits_;
};

using HexadecimalDeathTest = HexadecimalTest;

TEST_F(HexadecimalTest, EncodeAndDecode) {
  HexadecimalEncode(bytes_.get(), digits_.get());
  EXPECT_EQ(uppercase_digits_, digits_);
  UniqueBytes bytes(kBytes);
  HexadecimalDecode(digits_.get(), bytes.get());
  EXPECT_EQ(bytes_, bytes);
}

TEST_F(HexadecimalTest, InPlace) {
  auto buffer = std::make_unique<uint8_t[]>(kDigits);
  std::memcpy(&buffer[1], bytes_.data.get(), kBytes);
  HexadecimalEncode({&buffer[1], kBytes}, {&buffer[0], kDigits});
  EXPECT_EQ(uppercase_digits_, Bytes(&buffer[0], kDigits));
  std::memcpy(&buffer[0], bytes_.data.get(), kBytes);
  HexadecimalEncode({&buffer[0], kBytes}, {&buffer[0], kDigits});
  EXPECT_EQ(uppercase_digits_, Bytes(&buffer[0], kDigits));
  HexadecimalDecode({&buffer[0], kDigits}, {&buffer[1], kBytes});
  EXPECT_EQ(bytes_, Bytes(&buffer[1], kBytes));
  std::memcpy(&buffer[0], uppercase_digits_.data.get(), kDigits);
  HexadecimalDecode({&buffer[0], kDigits}, {&buffer[0], kBytes});
  EXPECT_EQ(bytes_, Bytes(&buffer[0], kBytes));
}

TEST_F(HexadecimalTest, LargeOutput) {
  int64_t const digits_size = kDigits + 42;
  UniqueBytes digits(digits_size);
  std::memset(digits.data.get(), 'X', digits_size);
  HexadecimalEncode(bytes_.get(), digits.get());
  EXPECT_EQ(uppercase_digits_, Bytes(digits.data.get(), kDigits));
  EXPECT_THAT(std::vector<uint8_t>(&digits.data[kDigits],
                                   &digits.data[digits_size]),
              Each('X'));
  int64_t const bytes_size = kBytes + 42;
  UniqueBytes bytes(bytes_size);
  std::memset(bytes.data.get(), 'Y', bytes_size);
  HexadecimalDecode(uppercase_digits_.get(), bytes.get());
  EXPECT_EQ(bytes_, Bytes(bytes.data.get(), kBytes));
  EXPECT_THAT(std::vector<uint8_t>(&bytes.data[kBytes],
                                   &bytes.data[bytes_size]),
              Each('Y'));
}

TEST_F(HexadecimalTest, Adjacent) {
  auto buffer = std::make_unique<uint8_t[]>(kDigits + kBytes);
  std::memcpy(&buffer[0], bytes_.data.get(), kBytes);
  HexadecimalEncode({&buffer[0], kBytes}, {&buffer[kBytes], kDigits});
  EXPECT_EQ(uppercase_digits_, Bytes(&buffer[kBytes], kDigits));
  std::memcpy(&buffer[0], uppercase_digits_.data.get(), kDigits);
  HexadecimalDecode({&buffer[0], kDigits}, {&buffer[kDigits], kBytes});
  EXPECT_EQ(bytes_, Bytes(&buffer[kDigits], kBytes));
  HexadecimalDecode({&buffer[0], kDigits + 1}, {&buffer[kDigits], kBytes});
  EXPECT_EQ(bytes_, Bytes(&buffer[kDigits], kBytes));
}

TEST_F(HexadecimalDeathTest, Overlap) {
  auto buffer = std::make_unique<uint8_t[]>(kDigits + kBytes - 1);
  EXPECT_DEATH({
    HexadecimalEncode({&buffer[kDigits - 1], kBytes}, {&buffer[0], kDigits});
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode({&buffer[0], kDigits}, {&buffer[kDigits - 1], kBytes});
  }, "bad overlap");
  buffer = std::make_unique<uint8_t[]>(kDigits);
  EXPECT_DEATH({
    HexadecimalEncode({&buffer[2], kBytes}, {&buffer[0], kDigits});
  }, "bad overlap");
  EXPECT_DEATH({
    HexadecimalDecode({&buffer[0], kDigits}, {&buffer[2], kBytes});
  }, "bad overlap");
}

TEST_F(HexadecimalDeathTest, Size) {
  std::vector<uint8_t> bytes(kBytes);
  std::vector<uint8_t> digits(kDigits);
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
  std::vector<uint8_t> bytes(kBytes);
  HexadecimalDecode(lowercase_digits_.get(), {bytes.data(), bytes.size()});
  EXPECT_EQ(bytes_, Bytes(bytes.data(), bytes.size()));
  HexadecimalDecode(uppercase_digits_.get(), {bytes.data(), bytes.size()});
  EXPECT_EQ(bytes_, Bytes(bytes.data(), bytes.size()));
}

TEST_F(HexadecimalTest, Invalid) {
  std::vector<uint8_t> bytes(1);
  std::vector<uint8_t> digits = {'a', 'b', 'c'};
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
