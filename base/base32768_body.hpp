
#pragma once

#include "base/base32768.hpp"

#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace base {

namespace internal_base32768 {

// Ceiling log2 of n.  8 -> 3, 7 -> 2.
constexpr int CeilingLog2(int const n) {
  return n == 1 ? 0 : CeilingLog2(n >> 1) + 1;
}

// A repertoire is used to encode or decode integers into Basic Multilingual
// Plane code points.
class Repertoire {
 public:
  virtual char16_t Encode(std::uint16_t k) const = 0;
  virtual std::uint16_t Decode(char16_t code_point) const = 0;
};

// A caching repertoire builds at constructions caches for fast encoding and
// decoding.
template<std::int64_t block_size, std::int64_t block_count>
class CachingRepertoire : public Repertoire {
 public:
  constexpr std::int64_t EncodingBits() const;

  // Returns true if this repertoire is capable of encoding the given code
  // point.
  bool CanEncode(char16_t code_point) const;

  char16_t Encode(std::uint16_t k) const override;
  std::uint16_t Decode(char16_t code_point) const override;

 private:
  template<std::int64_t block_count_plus_1>
  constexpr CachingRepertoire(
      char16_t const (&blocks)[block_count_plus_1]);

  char16_t const* const blocks_;
  std::int64_t const encoding_bits_;

  // These arrays are sparse: not all entries are filled with useful data.  The
  // caller must encode values which are within [0, block_count * block_size[,
  // and must decode characters which lie within the blocks of this repertoire.
  std::array<char16_t, block_count * block_size> encoding_cache_;
  // Using a C array because MSVC gets confused with an std::array.
  std::conditional_t<(CeilingLog2(block_size * block_count) > 8),
                     std::uint16_t,
                     std::uint8_t>
      decoding_cache_[std::numeric_limits<char16_t>::max()];

  template<std::int64_t s, std::int64_t c>
  friend constexpr CachingRepertoire<s, c - 1> MakeRepertoire(
    char16_t const (&blocks)[c]);
};

template<std::int64_t block_size, std::int64_t block_count_plus_1>
constexpr CachingRepertoire<block_size, block_count_plus_1 - 1> MakeRepertoire(
    char16_t const (&blocks)[block_count_plus_1]);

template<std::int64_t block_size, std::int64_t block_count>
template<std::int64_t block_count_plus_1>
constexpr CachingRepertoire<block_size, block_count>::CachingRepertoire(
    char16_t const (&blocks)[block_count_plus_1])
    : blocks_(blocks),
      encoding_bits_(CeilingLog2(block_size * block_count)),
      decoding_cache_() {
  // Don't do pointer arithmetic in this constructor, it confuses MSVC.
  static_assert(block_count_plus_1 == block_count + 1,
                "Incorrect literal size");

  // Check null-termination.
  assert(blocks_[block_count] == 0);
  // Check ordering and lack of overlap.
  for (int block = 0; block < block_count - 1; ++block) {
    assert(blocks_[block] < blocks_[block + 1]);
    assert(blocks_[block + 1] - blocks_[block] >= block_size);
  }

  for (int block = 0; block < block_count; ++block) {
    char16_t const block_start_code_point = blocks_[block];
    int const block_start_k = block_size * block;
    for (int offset = 0; offset < block_size; ++offset) {
      char16_t const code_point = block_start_code_point + offset;
      int const k = block_start_k + offset;
      encoding_cache_[k] = code_point;
      decoding_cache_[code_point] = k;
    }
  }
}

template<std::int64_t block_size, std::int64_t block_count>
constexpr std::int64_t
CachingRepertoire<block_size, block_count>::EncodingBits() const {
  return encoding_bits_;
}

template<std::int64_t block_size, std::int64_t block_count>
bool CachingRepertoire<block_size, block_count>::CanEncode(
    char16_t const code_point) const {
  char16_t const truncated_code_point = code_point & ~(block_size - 1);
  // Linear search because small: we only call this function for the 7-bit
  // encoding.
  for (int block = 0; block < block_count; ++block) {
    if (blocks_[block] == truncated_code_point) {
      return true;
    }
  }
  return false;
}

template<std::int64_t block_size, std::int64_t block_count>
char16_t CachingRepertoire<block_size, block_count>::Encode(
    std::uint16_t const k) const {
  // Check that the integer to encode has the expected number of bits.
  DCHECK_EQ(0, k & ~((1 << EncodingBits()) - 1)) << std::hex << k;
  return encoding_cache_[k];
}

template<std::int64_t block_size, std::int64_t block_count>
std::uint16_t CachingRepertoire<block_size, block_count>::Decode(
    char16_t const code_point) const {
  return decoding_cache_[code_point];
}

template<std::int64_t block_size, std::int64_t block_count_plus_1>
constexpr CachingRepertoire<block_size, block_count_plus_1 - 1> MakeRepertoire(
    char16_t const (&blocks)[block_count_plus_1]) {
  return CachingRepertoire<block_size, block_count_plus_1 - 1>(blocks);
}

constexpr int block_size = 1 << 5;

constexpr auto fifteen_bits = MakeRepertoire<block_size>(
    u"ҠԀڀڠݠހ߀ကႠᄀᄠᅀᆀᇠሀሠበዠጠᎠᏀᐠᑀᑠᒀᒠᓀᓠᔀᔠᕀᕠ"
    u"ᖀᖠᗀᗠᘀᘠᙀᚠᛀកᠠᡀᣀᦀ᧠ᨠᯀᰀᴀ⇠⋀⍀⍠⎀⎠⏀␀─┠╀╠▀"
    u"■◀◠☀☠♀♠⚀⚠⛀⛠✀✠❀➀➠⠀⠠⡀⡠⢀⢠⣀⣠⤀⤠⥀⥠⦠⨠⩀⪀"
    u"⪠⫠⬀⬠⭀ⰀⲀⲠⳀⴀⵀ⺠⻀㇀㐀㐠㑀㑠㒀㒠㓀㓠㔀㔠㕀㕠㖀㖠㗀㗠㘀㘠"
    u"㙀㙠㚀㚠㛀㛠㜀㜠㝀㝠㞀㞠㟀㟠㠀㠠㡀㡠㢀㢠㣀㣠㤀㤠㥀㥠㦀㦠㧀㧠㨀㨠"
    u"㩀㩠㪀㪠㫀㫠㬀㬠㭀㭠㮀㮠㯀㯠㰀㰠㱀㱠㲀㲠㳀㳠㴀㴠㵀㵠㶀㶠㷀㷠㸀㸠"
    u"㹀㹠㺀㺠㻀㻠㼀㼠㽀㽠㾀㾠㿀㿠䀀䀠䁀䁠䂀䂠䃀䃠䄀䄠䅀䅠䆀䆠䇀䇠䈀䈠"
    u"䉀䉠䊀䊠䋀䋠䌀䌠䍀䍠䎀䎠䏀䏠䐀䐠䑀䑠䒀䒠䓀䓠䔀䔠䕀䕠䖀䖠䗀䗠䘀䘠"
    u"䙀䙠䚀䚠䛀䛠䜀䜠䝀䝠䞀䞠䟀䟠䠀䠠䡀䡠䢀䢠䣀䣠䤀䤠䥀䥠䦀䦠䧀䧠䨀䨠"
    u"䩀䩠䪀䪠䫀䫠䬀䬠䭀䭠䮀䮠䯀䯠䰀䰠䱀䱠䲀䲠䳀䳠䴀䴠䵀䵠䶀䷀䷠一丠乀"
    u"习亀亠什仠伀传佀你侀侠俀俠倀倠偀偠傀傠僀僠儀儠兀兠冀冠净几刀删剀"
    u"剠劀加勀勠匀匠區占厀厠叀叠吀吠呀呠咀咠哀哠唀唠啀啠喀喠嗀嗠嘀嘠噀"
    u"噠嚀嚠囀因圀圠址坠垀垠埀埠堀堠塀塠墀墠壀壠夀夠奀奠妀妠姀姠娀娠婀"
    u"婠媀媠嫀嫠嬀嬠孀孠宀宠寀寠尀尠局屠岀岠峀峠崀崠嵀嵠嶀嶠巀巠帀帠幀"
    u"幠庀庠廀廠开张彀彠往徠忀忠怀怠恀恠悀悠惀惠愀愠慀慠憀憠懀懠戀戠所"
    u"扠技抠拀拠挀挠捀捠掀掠揀揠搀搠摀摠撀撠擀擠攀攠敀敠斀斠旀无昀映晀"
    u"晠暀暠曀曠最朠杀杠枀枠柀柠栀栠桀桠梀梠检棠椀椠楀楠榀榠槀槠樀樠橀"
    u"橠檀檠櫀櫠欀欠歀歠殀殠毀毠氀氠汀池沀沠泀泠洀洠浀浠涀涠淀淠渀渠湀"
    u"湠満溠滀滠漀漠潀潠澀澠激濠瀀瀠灀灠炀炠烀烠焀焠煀煠熀熠燀燠爀爠牀"
    u"牠犀犠狀狠猀猠獀獠玀玠珀珠琀琠瑀瑠璀璠瓀瓠甀甠畀畠疀疠痀痠瘀瘠癀"
    u"癠皀皠盀盠眀眠着睠瞀瞠矀矠砀砠础硠碀碠磀磠礀礠祀祠禀禠秀秠稀稠穀"
    u"穠窀窠竀章笀笠筀筠简箠節篠簀簠籀籠粀粠糀糠紀素絀絠綀綠緀締縀縠繀"
    u"繠纀纠绀绠缀缠罀罠羀羠翀翠耀耠聀聠肀肠胀胠脀脠腀腠膀膠臀臠舀舠艀"
    u"艠芀芠苀苠茀茠荀荠莀莠菀菠萀萠葀葠蒀蒠蓀蓠蔀蔠蕀蕠薀薠藀藠蘀蘠虀"
    u"虠蚀蚠蛀蛠蜀蜠蝀蝠螀螠蟀蟠蠀蠠血衠袀袠裀裠褀褠襀襠覀覠觀觠言訠詀"
    u"詠誀誠諀諠謀謠譀譠讀讠诀诠谀谠豀豠貀負賀賠贀贠赀赠趀趠跀跠踀踠蹀"
    u"蹠躀躠軀軠輀輠轀轠辀辠迀迠退造遀遠邀邠郀郠鄀鄠酀酠醀醠釀釠鈀鈠鉀"
    u"鉠銀銠鋀鋠錀錠鍀鍠鎀鎠鏀鏠鐀鐠鑀鑠钀钠铀铠销锠镀镠門閠闀闠阀阠陀"
    u"陠隀隠雀雠需霠靀靠鞀鞠韀韠頀頠顀顠颀颠飀飠餀餠饀饠馀馠駀駠騀騠驀"
    u"驠骀骠髀髠鬀鬠魀魠鮀鮠鯀鯠鰀鰠鱀鱠鲀鲠鳀鳠鴀鴠鵀鵠鶀鶠鷀鷠鸀鸠鹀"
    u"鹠麀麠黀黠鼀鼠齀齠龀龠ꀀꀠꁀꁠꂀꂠꃀꃠꄀꄠꅀꅠꆀꆠꇀꇠꈀꈠꉀꉠꊀ"
    u"ꊠꋀꋠꌀꌠꍀꍠꎀꎠꏀꏠꐀꐠꑀꑠ꒠ꔀꔠꕀꕠꖀꖠꗀꗠꙀꚠꛀ꜀꜠ꝀꞀꡀ");
constexpr auto seven_bits = MakeRepertoire<block_size>(u"ƀɀɠʀ");

constexpr std::int64_t bits_per_byte = 8;
constexpr std::int64_t bits_per_code_point = fifteen_bits.EncodingBits();
constexpr std::int64_t bits_per_final_code_point = seven_bits.EncodingBits();
constexpr std::int64_t bytes_per_code_point =
    (bits_per_code_point + 2 * bits_per_byte - 2) / bits_per_byte;
static_assert(bytes_per_code_point == 3,
              "End of input padding below won't be correct");

void Base32768Encode(Array<std::uint8_t const> input,
                     Array<char16_t> output) {
  CHECK_NOTNULL(input.data);
  CHECK(input.size == 0 || output.data != nullptr);

  std::uint8_t const* const input_end = input.data + input.size;
  std::int64_t input_bit_index = 0;
  while (input.data < input_end) {
    std::int32_t data;

    // Prepare for normal encoding.
    std::int32_t shift = bytes_per_code_point * bits_per_byte -
                         bits_per_code_point - input_bit_index;
    std::int32_t mask = ((1 << bits_per_code_point) - 1) << shift;
    Repertoire const* repertoire = &fifteen_bits;

    if (input_end - input.data >= bytes_per_code_point) {
      // Extract three bytes.
      data = input.data[0] << 2 * bits_per_byte |
             input.data[1] << bits_per_byte |
             input.data[2];
    } else if (input_end - input.data == 2) {
      // Extract the last two bytes and pad with 1s.
      data = input.data[0] << 2 * bits_per_byte |
             input.data[1] << bits_per_byte |
             ((1 << bits_per_byte) - 1);
    } else if (input_end - input.data == 1) {
      // Extract the last byte and pad with 1s.
      data = input.data[0] << 2 * bits_per_byte |
             ((1 << 2 * bits_per_byte) - 1);
      if (input_bit_index > 0) {
        // Switch to special encoding.
        shift = bytes_per_code_point * bits_per_byte -
                bits_per_final_code_point - input_bit_index;
        mask = ((1 << bits_per_final_code_point) - 1) << shift;
        repertoire = &seven_bits;
      }
    }
    std::int32_t code_point = (data & mask) >> shift;
    CHECK_LE(0, code_point);
    CHECK_LT(code_point, 1 << bits_per_code_point);
    output.data[0] = repertoire->Encode(code_point);

    // The following computation may cause |input.data| to overshoot the end if
    // using the special encoding at the end.  This is safe as soon as the loop
    // condition uses <.
    input_bit_index += bits_per_code_point;
    input.data += input_bit_index / bits_per_byte;
    input_bit_index %= bits_per_byte;
    ++output.data;
  }
}
UniqueArray<char16_t> Base32768Encode(Array<std::uint8_t const> input,
                                      bool const null_terminated) {
  base::UniqueArray<char16_t> output(Base32768EncodedLength(input) +
                                     (null_terminated ? 1 : 0));
  if (output.size > 0) {
    base::Base32768Encode(input, output.get());
  }
  if (null_terminated) {
    output.data[output.size - 1] = 0;
  }
  return output;
}

std::int64_t Base32768EncodedLength(Array<std::uint8_t const> const input) {
  return (input.size * bits_per_byte + bits_per_code_point - 1) /
         bits_per_code_point;
}

void Base32768Decode(Array<char16_t const> input, Array<std::uint8_t> output) {
  CHECK_NOTNULL(input.data);
  CHECK(input.size == 0 || output.data != nullptr);

  char16_t const* const input_end = input.data + input.size;
  std::uint8_t const* const output_end = output.data + output.size;
  std::int64_t output_bit_index = 0;
  while (input.data < input_end) {
    bool const at_end = input_end - input.data == 1;
    std::int32_t data;
    std::int32_t shift = bytes_per_code_point * bits_per_byte -
                         bits_per_code_point - output_bit_index;
    Repertoire const* repertoire = &fifteen_bits;
    if (at_end && seven_bits.CanEncode(input.data[0])) {
      shift = bytes_per_code_point * bits_per_byte - bits_per_final_code_point -
              output_bit_index;
      repertoire = &seven_bits;
    }

    // Align |data| on the output bit index.
    data = repertoire->Decode(input.data[0]);
    data <<= shift;

    // Fill the output with the parts of the code point belonging to each byte.
    if (output_bit_index == 0) {
      output.data[0] = (data >> (2 * bits_per_byte));
    } else {
      output.data[0] |= (data >> (2 * bits_per_byte));
    }
    if (shift < 2 * bits_per_byte && output_end - output.data > 1) {
      output.data[1] = (data >> bits_per_byte) & ((1 << bits_per_byte) - 1);
      if (shift < bits_per_byte && output_end - output.data > 2) {
        output.data[2] = data & ((1 << bits_per_byte) - 1);
      }
    }

    output_bit_index += bits_per_code_point;
    output.data += output_bit_index / bits_per_byte;
    output_bit_index %= bits_per_byte;
    ++input.data;
  }
}

UniqueArray<std::uint8_t> Base32768Decode(Array<char16_t const> input) {
  UniqueArray<std::uint8_t> output(Base32768DecodedLength(input));
  if (output.size > 0) {
    Base32768Decode(input, output.get());
  }
  return output;
}

std::int64_t Base32768DecodedLength(Array<char16_t const> const input) {
  // In order to decide how many bytes the input will decode to, we need to
  // figure out if the last code point encodes 7 bits or 15 bits.
  std::int64_t encoded_bits;
  if (input.size > 0 && seven_bits.CanEncode(input.data[input.size - 1])) {
    encoded_bits =
        (input.size - 1) * bits_per_code_point + bits_per_final_code_point;
  } else {
    encoded_bits = input.size * bits_per_code_point;
  }
  // Either we have a multiple of 15 bits, in which case the division is exact;
  // or there is padding, and truncation has the right effect.
  return encoded_bits / bits_per_byte;
}

}  // namespace internal_base32768
}  // namespace base
}  // namespace principia
