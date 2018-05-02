
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
    char16_t const* const output32768 =
        reinterpret_cast<char16_t*>(output.get().data);
    EXPECT_EQ(0,
              std::char_traits<char16_t>::compare(
                  output32768, base32768, base32768_length));
    LOG(ERROR) << base32768_length;
    for (int i = 0; i < base32768_length; ++i) {
      EXPECT_EQ(output32768[i], base32768[i])
          << "index=" << i << " actual=" << std::hex << output32768[i]
          << " expected=" << base32768[i];
    }
  }
};

using Base32768DeathTest = Base32768Test;

TEST_F(Base32768Test, EncodeCaseDemo) {
  char const binary[] = "\xd4\x1d\x8c\xd9\x8f\x00\xb2\x04"
                        "\xe9\x80\x09\x98\xec\xf8\x42";
  char16_t const base32768[] = u"遮視塀⤠䶌Ԇ堹麢";
  CheckEncoding(binary, /*binary_length=*/15, base32768);
}

TEST_F(Base32768Test, EncodeCaseEmpty) {
  char const binary[] = "";
  char16_t const base32768[] = u"";
  CheckEncoding(binary, /*binary_length=*/0, base32768);
}

TEST_F(Base32768Test, EncodeEveryByte) {
  char const binary[] = "\x00\x01\x02\x03\x04\x05\x06\x07"
                        "\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f"
                        "\x10\x11\x12\x13\x14\x15\x16\x17"
                        "\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f"
                        "\x20\x21\x22\x23\x24\x25\x26\x27"
                        "\x28\x29\x2a\x2b\x2c\x2d\x2e\x2f"
                        "\x30\x31\x32\x33\x34\x35\x36\x37"
                        "\x38\x39\x3a\x3b\x3c\x3d\x3e\x3f"
                        "\x40\x41\x42\x43\x44\x45\x46\x47"
                        "\x48\x49\x4a\x4b\x4c\x4d\x4e\x4f"
                        "\x50\x51\x52\x53\x54\x55\x56\x57"
                        "\x58\x59\x5a\x5b\x5c\x5d\x5e\x5f"
                        "\x60\x61\x62\x63\x64\x65\x66\x67"
                        "\x68\x69\x6a\x6b\x6c\x6d\x6e\x6f"
                        "\x70\x71\x72\x73\x74\x75\x76\x77"
                        "\x78\x79\x7a\x7b\x7c\x7d\x7e\x7f"
                        "\x80\x81\x82\x83\x84\x85\x86\x87"
                        "\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f"
                        "\x90\x91\x92\x93\x94\x95\x96\x97"
                        "\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f"
                        "\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7"
                        "\xa8\xa9\xaa\xab\xac\xad\xae\xaf"
                        "\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7"
                        "\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf"
                        "\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7"
                        "\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf"
                        "\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7"
                        "\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf"
                        "\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7"
                        "\xe8\xe9\xea\xeb\xec\xed\xee\xef"
                        "\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7"
                        "\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff";
  char16_t const base32768[] = u"Ҡ曠蛠盀庠䩨㱘Ⳏ"
                               u"┨ᗄ棂枱團蛄媖䉝"
                               u"㕏湨䪄墢侑䋠碴册"
                               u"㳗⩌ᰆ䥳䟩缼雒悛"
                               u"䑞痯蹨㩤䁢㭙㓐澪"
                               u"䯦㦓灊ᛕ㣚瞵匎纹"
                               u"厍絷别鱦⤓⻑焬跈"
                               u"嬕䄛㇍赗ᔋ瀭轊鳗"
                               u"抜蓾闯繇ꊤᯊ⎩ᣆ"
                               u"樤䢢矑漸髜梦䭧㫕"
                               u"熫貆妳怩鍔ꖂ榥䧤"
                               u"礳偊㭵儚词愞蟃夓"
                               u"肺鐍鵷䇫葅鵚꜡栢"
                               u"衂埑罙⭜精妗䏟眱"
                               u"迉鮕愺ꐭ甶闓戝虀"
                               u"靑彙䋼铞涯刏耻镏"
                               u"默ꍜꖞ藏昧蹋鹙꒾"
                               u"ꡟ";
  CheckEncoding(binary, /*binary_length=*/256, base32768);
}

}  // namespace base
}  // namespace principia
