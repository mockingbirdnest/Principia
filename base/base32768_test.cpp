
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
    Bytes input({reinterpret_cast<std::uint8_t const*>(binary), binary_length});
    UniqueBytes output(base32768_length);
    Base32768Encode(input, output.get());
    EXPECT_EQ(0,
              std::char_traits<char16_t>::compare(
                  output.get().data, base32768, base32768_length));
  }
};

using Base32768DeathTest = Base32768Test;

TEST_F(Base32768Test, EncodeAndDecode) {
  char const binary[] = "ÔŒÙ ²é€	˜ìøB";
  char16_t const base32768[] = u"遮視塀⤠䶌Ԇ堹麢";
}

}  // namespace base
}  // namespace principia
