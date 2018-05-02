
#include "base/base32768.hpp"

#include <memory>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace base {

class Base32768Test : public testing::Test {
 protected:
  void CheckEncoding(char const* const binary,
                     int const binary_length,
                     const char16_t* const base32768) {
    int const base32768_length = std::char_traits<char16_t>::length(base32768);
    Array<std::uint8_t const> input(
        reinterpret_cast<std::uint8_t const*>(binary), binary_length);
    UniqueBytes output(sizeof(char16_t) * base32768_length);
    Base32768Encode(input, output.get());
    EXPECT_EQ(0,
              std::char_traits<char16_t>::compare(
                  reinterpret_cast<char16_t*>(output.get().data),
                  base32768,
                  base32768_length));
  }
};

using Base32768DeathTest = Base32768Test;

TEST_F(Base32768Test, EncodeAndDecode) {
  char const binary[] = "\xd4\x1d\x8c\xd9\x8f\x00\xb2\x04"
                        "\xe9\x80\x09\x98\xec\xf8\x42";
  char16_t const base32768[] = u"遮視塀⤠䶌Ԇ堹麢";
  CheckEncoding(binary, /*binary_length=*/15, base32768);
}

}  // namespace base
}  // namespace principia
