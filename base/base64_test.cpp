
#include "base/base64.hpp"

#include <string>

#include "gtest/gtest.h"

namespace principia {
namespace base {

class Base64Test : public ::testing::Test {
 protected:
  Base64Encoder</*null_terminated=*/false> encoder_;
};

TEST_F(Base64Test, WellKnownData) {
  // From https://en.wikipedia.org/wiki/Base64.
  std::string const decoded_string =
      "Man is distinguished, not only by his reason, but by this singular "
      "passion from other animals, which is a lust of the mind, that by a "
      "perseverance of delight in the continued and indefatigable generation "
      "of knowledge, exceeds the short vehemence of any carnal pleasure.";
  std::string const encoded_string =
      "TWFuIGlzIGRpc3Rpbmd1aXNoZWQsIG5vdCBvbmx5IGJ5IGhpcyByZWFzb24sIGJ1dCBieSB0"
      "aGlzIHNpbmd1bGFyIHBhc3Npb24gZnJvbSBvdGhlciBhbmltYWxzLCB3aGljaCBpcyBhIGx1"
      "c3Qgb2YgdGhlIG1pbmQsIHRoYXQgYnkgYSBwZXJzZXZlcmFuY2Ugb2YgZGVsaWdodCBpbiB0"
      "aGUgY29udGludWVkIGFuZCBpbmRlZmF0aWdhYmxlIGdlbmVyYXRpb24gb2Yga25vd2xlZGdl"
      "LCBleGNlZWRzIHRoZSBzaG9ydCB2ZWhlbWVuY2Ugb2YgYW55IGNhcm5hbCBwbGVhc3VyZS4";

  {
    Array<std::uint8_t const> decoded_array(
        reinterpret_cast<std::uint8_t const*>(decoded_string.c_str()),
        decoded_string.size());
    auto const encoded_array = encoder_.Encode(decoded_array);
    EXPECT_EQ(encoded_string.size(), encoder_.EncodedLength(decoded_array));
    for (std::int64_t i = 0; i < encoded_array.size; ++i) {
      EXPECT_EQ(encoded_string[i], encoded_array.data[i]) << i;
    }
  }
  {
    Array<char const> encoded_array(encoded_string.c_str(),
                                    encoded_string.size());
    auto const decoded_array = encoder_.Decode(encoded_array);
    EXPECT_EQ(decoded_string.size(), encoder_.DecodedLength(encoded_array));
    for (std::int64_t i = 0; i < decoded_array.size; ++i) {
      EXPECT_EQ(decoded_string[i],
                static_cast<char>(decoded_array.data[i])) << i;
    }
  }
}

}  // namespace base
}  // namespace principia
