
#include "base/base32768.hpp"

#include <memory>
#include <random>
#include <string>
#include <vector>

#include "base/array.hpp"
#include "base/macros.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// Clang doesn't have a correct |std::array| yet, and we don't actually use this
// code, so let's get rid of the entire test.
#if PRINCIPIA_COMPILER_MSVC

MATCHER_P(EqualsBytes, expected, "") {
  auto const actual = arg;
  if (actual.size != expected.size) {
    *result_listener << "\nDifferent length: " << actual.size << " vs. "
                     << expected.size;
    return false;
  }
  for (int i = 0; i < actual.size; ++i) {
    if (actual.data[i] != expected.data[i]) {
      *result_listener << "\nBytes differ at index " << i << "/" << actual.size
                       << ": " << std::hex << static_cast<int>(actual.data[i])
                       << " vs. " << static_cast<int>(expected.data[i]);
      return false;
    }
  }
  return true;
}

namespace principia {
namespace base {

void PrintTo(Array<std::uint8_t> bytes, std::ostream* os) {
  *os << std::hex;
  for (int i = 0; i < bytes.size; ++i) {
    *os << static_cast<int>(bytes.data[i]) << " ";
  }
  *os << "\n";
}

class Base32768Test : public testing::Test {
 protected:
  void CheckEncoding(Array<std::uint8_t const> const binary,
                     Array<char16_t const> const base32768) {
    CHECK_EQ(base32768.size, encoder_.EncodedLength(binary));
    UniqueArray<char16_t> output(base32768.size);
    encoder_.Encode(binary, output.get());
    EXPECT_EQ(0,
              std::char_traits<char16_t>::compare(
                  output.data.get(), base32768.data, base32768.size));
    for (int i = 0; i < base32768.size; ++i) {
      EXPECT_EQ(output.data[i], base32768.data[i])
          << "index=" << i << " actual=" << std::hex << output.data[i]
          << " expected=" << base32768.data[i];
    }
  }

  void CheckDecoding(Array<std::uint8_t const> const binary,
                     Array<char16_t const> const base32768) {
    CHECK_EQ(binary.size, encoder_.DecodedLength(base32768));
    UniqueArray<std::uint8_t> output(binary.size);
    encoder_.Decode(base32768, output.get());
    EXPECT_EQ(0, std::memcmp(output.data.get(), binary.data, binary.size));
    for (int i = 0; i < binary.size; ++i) {
      EXPECT_EQ(output.data[i], binary.data[i]) << "index=" << i;
    }
  }

  Base32768Encoder</*null_terminated=*/false> encoder_;
};

using Base32768DeathTest = Base32768Test;

// No tests because of a bug in 15.8 preview 3.
#if !PRINCIPIA_COMPILER_MSVC || \
    !(_MSC_FULL_VER == 191526608 || \
      _MSC_FULL_VER == 191526731 || \
      _MSC_FULL_VER == 191627024)
TEST_F(Base32768Test, EncodeMultipleOf15Bits) {
  // First 15 bytes of the MD5 of the empty string.
  Array<std::uint8_t const> const binary("\xd4\x1d\x8c\xd9\x8f\x00\xb2\x04"
                                         "\xe9\x80\x09\x98\xec\xf8\x42");
  Array<char16_t const> const base32768(u"遮視塀⤠䶌Ԇ堹麢");
  CheckEncoding(binary, base32768);
  CheckDecoding(binary, base32768);
}

TEST_F(Base32768Test, EncodeEmpty) {
  Array<std::uint8_t const> const binary("");
  Array<char16_t const> const base32768(u"");
  CheckEncoding(binary, base32768);
  CheckDecoding(binary, base32768);
}

TEST_F(Base32768Test, EncodeEveryByte) {
  std::array<std::uint8_t, 256> binary;
  for (int c = 0; c < binary.size(); ++c) {
    binary[c] = c;
  }
  Array<char16_t const> const base32768(u"Ҡ曠蛠盀庠䩨㱘Ⳏ"
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
                                        u"ꡟ");
  CheckEncoding(binary, base32768);
  CheckDecoding(binary, base32768);
}

TEST_F(Base32768Test, EncodeHatetrisWrRle) {
  Array<std::uint8_t const> const binary("\xc0\x03\x91\xc0\x03\x8d\xc0\x02"
                                         "\x80\x00\x8b\xc0\x01\x90\xc0\x01"
                                         "\x8e\xc0\x90\x42\x91\xc0\x43\x8e"
                                         "\x41\x90\xc0\x42\x8e\x40\x90\xc0"
                                         "\x41\x9f\xc0\x40\x8e\x41\xc0\x42"
                                         "\x8d\x00\x90\xc0\x8e\xc0\x41\x8c"
                                         "\xc0\x42\x97\x01\x8b\xc0\x03\x8b"
                                         "\x40\x80\xc0\x41\x83\xc0\x44\x8d"
                                         "\x41\x00\x8b\x40\x8a\xc0\x43\x89"
                                         "\xc0\x42\x8a\x02\xc0\x00\x8c\x40"
                                         "\x80\x41\x80\xc2\x03\x8d\x40\x80"
                                         "\xc0\x80\xc0\x80\xc2\x03\x8e\xc0"
                                         "\x81\xc0\x03\x8e\x02\x8e\x03\x8e"
                                         "\x42\x00\x8d\x00\x40\x8d\x01\x8d"
                                         "\x03\x8d\xc0\x03\x8b\xc0\x04\x89"
                                         "\xc0\x03\x87\xc0\x01\x8b\xc0\x01"
                                         "\x89\xc0\x8b\x42\x8c\x40\xc0\x42"
                                         "\x89\x41\x8c\xc0\x42\x88\x40\x8b"
                                         "\xc0\x41\x95\xc0\x40\x89\xc0\x45"
                                         "\x88\x00\x8b\xc0\x89\xc0\x41\x87"
                                         "\xc0\x42\x8d\x01\x86\xc0\x03\x86"
                                         "\x40\xc0\x80\xc0\x41\x83\xc0\x45"
                                         "\x88\x41\x00\x86\x40\x85\xc0\x43"
                                         "\x84\xc0\x42\x85\x02\xc2\x00\x89"
                                         "\x43\x80\x01\xc2\x01\x88\x40\x80"
                                         "\xc0\x80\xc0\x80\xc2\x03\x89\xc0"
                                         "\x81\x03\xc0\x00\x89\x01\x80\xc1"
                                         "\x88\x03\x89\x41\x88\x00\x40\x88"
                                         "\x01\x88\x03\x88\xc0\x03\x86\xc0"
                                         "\x03\x84\xc0\x01\x86\xc0\x01\x84"
                                         "\xc0\x86\x42\x87\xc0\x43\x84\x41"
                                         "\x86\xc0\x42\x84\xc0\x03\x83\x40"
                                         "\x86\xc0\x41\x8b\xc0\x40\x84\xc0"
                                         "\x43\x83\x00\x86\xc0\x84\x01\x86"
                                         "\xc0\x00\x84\xc0\x43\x83\xc0\x41"
                                         "\x83\x02\x86\x03\x86\xc0\x43\x83"
                                         "\xc0\x41\x83\xc0\x83\xc0\x01\x85"
                                         "\xc0\x01\x83\x01\xc0\x01\x84\xc0"
                                         "\x03\x82");
  Array<char16_t const> const base32768(u"虡諐ڱ艠㩀ዯң㜀߇"
                                        u"嚄◲䉄䊲ᴃݥ㒀滀"
                                        u"㚧鹨⚄㑂⠴ᄁ曮蚀"
                                        u"覐◲陸葠㑭Ⴁ暡桀"
                                        u"㝣亀♔ᗖႮ㧀棪ᅀ"
                                        u"ԃ▰ᖘᯐ㑵Ⴁ曠蚠"
                                        u"因頸▼Ҽ幪䉇㒂چ"
                                        u"晰㟠㼐䊮Ү㷀ᘉ虡衐"
                                        u"ԑ扠⬎ይᛘ朠䞄盃"
                                        u"㹈令ᗞႦ几曩蚂衠"
                                        u"㞸☼ቬ䕁ᛚᆆ虡蟰"
                                        u"㹐⪤⪾Ⴖ㛂ݦ䚂陰"
                                        u"雸ᖨ乶▢㣇ҡ蝠衰"
                                        u"㙘▬ᖆ▮㧁Ⴃ虠䢀"
                                        u"噸㻀䊊ᰀݡ□橡袐"
                                        u"ڰ鉠䉦Ҧ⻀ᆄ蚣㛡"
                                        u"鹨庤⫶Ⴊ✠ᔃ䚃噰"
                                        u"埘ᖈ䱂㑌Ⴍ曤߃噠"
                                        u"㛘ᖸ䑂ᯌᛌᔆ蚁蝐"
                                        u"◐扨䑀ᯗңᒁ虠螐ڰ"
                                        u"ɏ");
  CheckEncoding(binary, base32768);
  CheckDecoding(binary, base32768);
}

TEST_F(Base32768Test, EncodeHatetrisWr) {
  Array<std::uint8_t const> const binary("\xc0\x2a\xaa\xaa\xaa\xab\x00\xaa"
                                         "\xaa\xaa\xac\x08\xaa\xaa\xaa\xc2"
                                         "\xaa\xaa\xaa\xaa\xc2\xaa\xaa\xaa"
                                         "\xae\xaa\xaa\xaa\xaa\x56\xaa\xaa"
                                         "\xaa\xaa\xb5\x5a\xaa\xaa\xaa\x96"
                                         "\xaa\xaa\xaa\xaa\xd5\xaa\xaa\xaa"
                                         "\xa9\xaa\xaa\xaa\xaa\xb5\xaa\xaa"
                                         "\xaa\xaa\xaa\xaa\xaa\xaa\xda\xaa"
                                         "\xaa\xaa\x97\x56\xaa\xaa\xaa\x8a"
                                         "\xaa\xaa\xaa\xab\xaa\xaa\xaa\xab"
                                         "\x5a\xaa\xaa\xab\x56\xaa\xaa\xaa"
                                         "\xaa\xaa\xa8\x2a\xaa\xaa\xb0\x0a"
                                         "\xaa\xaa\xa6\xd6\xab\x55\x6a\xaa"
                                         "\xaa\xa9\x4a\xaa\xaa\xa6\xaa\xaa"
                                         "\xad\x56\xaa\xaa\xb5\x6a\xaa\xaa"
                                         "\x03\x2a\xaa\xaa\xa6\x5b\xf0\x0a"
                                         "\xaa\xaa\xaa\x6e\xef\xc0\x2a\xaa"
                                         "\xaa\xaa\xeb\x00\xaa\xaa\xaa\xa8"
                                         "\x0a\xaa\xaa\xaa\x80\x2a\xaa\xaa"
                                         "\xaa\x54\xaa\xaa\xaa\xa1\xaa\xaa"
                                         "\xaa\xa0\xaa\xaa\xaa\xa0\x0a\xaa"
                                         "\xaa\xaa\xc0\x2a\xaa\xaa\xb0\x02"
                                         "\xaa\xaa\xb0\x0a\xaa\xac\x2a\xaa"
                                         "\xaa\xb0\xaa\xaa\xae\xaa\xaa\xa9"
                                         "\x5a\xaa\xaa\xa9\xd5\xaa\xaa\xa5"
                                         "\xaa\xaa\xaa\xb5\x6a\xaa\xa6\xaa"
                                         "\xaa\xab\x5a\xaa\xaa\xaa\xaa\xaa"
                                         "\xda\xaa\xaa\xd5\x56\xaa\xaa\x2a"
                                         "\xaa\xaa\xba\xaa\xaa\xd6\xaa\xab"
                                         "\x56\xaa\xaa\xaa\x82\xaa\xac\x02"
                                         "\xaa\xa7\xb5\xaa\xd5\x56\xaa\xaa"
                                         "\x52\xaa\xa6\xaa\xb5\x5a\xab\x56"
                                         "\xaa\x80\xfc\xaa\xaa\xa5\x58\x3f"
                                         "\x0a\xaa\xa9\xbb\xbf\x00\xaa\xaa"
                                         "\xae\x80\x32\xaa\xaa\x82\xfa\xaa"
                                         "\xa8\x02\xaa\xaa\x96\xaa\xaa\x1a"
                                         "\xaa\xa8\x2a\xaa\xa0\x0a\xaa\xab"
                                         "\x00\xaa\xab\x00\xaa\xb0\xaa\xab"
                                         "\x0a\xab\xaa\xa9\x5a\xaa\xad\x56"
                                         "\xaa\x5a\xaa\xb5\x6a\xac\x02\xa9"
                                         "\xaa\xab\x5a\xaa\xaa\xad\xaa\xb5"
                                         "\x5a\xa2\xaa\xae\xaa\x0a\xaa\xb2"
                                         "\xaa\xd5\x6a\xb5\xaa\x02\xaa\xa0"
                                         "\x0a\xaa\xb5\x5a\xad\x6a\xba\xac"
                                         "\x2a\xab\x0a\xa0\xc2\xaa\xc0\x2a");
  Array<char16_t const> const base32768(u"虵儊箵噪箵儐㞕儊"
                                        u"螵儊箸儊箵愊箵傶"
                                        u"箵儊紋儊箴脊箵儵"
                                        u"箵儊宵儊簍儊箵儊"
                                        u"箵崊箵俕宵儊㮕儊"
                                        u"簵儊篋儊箺脊箵儊"
                                        u"穵儊籡儊箖脍儵儊"
                                        u"笅儊笵儊鄕儊鄵儊"
                                        u"ᆕ儊笫敠箵儉萿暊"
                                        u"箵儚虵儊箠儊箵Ԋ"
                                        u"箵僵㮕儊ⵕ儊枵儊"
                                        u"癥儊箸ᐪ箵晪箵噪"
                                        u"箶⢪箶⢪篕儊礕儊"
                                        u"筊鄊笕儊簋儊玵"
                                        u"儋厵儊箵吊箶箺箴"
                                        u"儊箽儊脵儕宵儊ᠵ"
                                        u"剢箳鏊脊鄊磵僺篊"
                                        u"脋况䙿㮕債桘儊秗"
                                        u"敠箵刀ᴕ儈㸕儀ᠵ"
                                        u"僶箵⍊筥儊ڕ儌ᄵ"
                                        u"兠箸儊螵愊焵儕厴"
                                        u"脊脕兠箕儕箵儖篊"
                                        u"脂箷僢箶儍况紈ᠵ"
                                        u"䙊箺紋厷儢箸僣"
                                        u"ᠵ暊");
  CheckEncoding(binary, base32768);
  CheckDecoding(binary, base32768);
}

TEST_F(Base32768Test, Random) {
  std::mt19937_64 random(42);
  std::uniform_int_distribution<std::uint64_t> length_distribution(100, 150);
  std::uniform_int_distribution<int> bytes_distribution(0, 256);
  for (int test = 0; test < 1000; ++test) {
    UniqueArray<std::uint8_t> binary1(length_distribution(random));
    for (int i = 0; i < binary1.size; ++i) {
      binary1.data[i] = bytes_distribution(random);
    }

    UniqueArray<char16_t> const base32768 =
        Base32768Encode(binary1.get(),
                        /*null_terminated=*/false);
    UniqueArray<std::uint8_t> binary2 = Base32768Decode(base32768.get());

    EXPECT_THAT(binary2.get(), EqualsBytes(binary1.get())) << "test: " << test;
  }
}
#endif

}  // namespace base
}  // namespace principia

#endif
