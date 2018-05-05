
#pragma once

#include "base/base32768.hpp"

#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace base {

namespace internal_base32768 {

constexpr std::int64_t bits_per_byte = 8;
constexpr std::int64_t bits_per_code_point = 15;
constexpr std::int64_t bits_per_final_code_point = 7;
constexpr std::int64_t bytes_per_code_point =
    (bits_per_code_point + 2 * bits_per_byte - 2) / bits_per_byte;
static_assert(bytes_per_code_point == 3,
              "End of input padding below won't be correct");

constexpr int block_size = 1 << 5;

class Repertoire {
 public:
  //constexpr Repertoire() = default;
  //virtual ~Repertoire() = default;

  virtual char16_t const& Encode(std::uint16_t const k) const = 0;
  virtual std::uint16_t const& Decode(char16_t const code_point) const = 0;
};

template<std::size_t block_size, std::size_t block_count>
class CachingRepertoire : public Repertoire {
 public:
  template<std::size_t block_count_plus_1>
  constexpr CachingRepertoire(
      char16_t const (&blocks_begin)[block_count_plus_1]);

  char16_t const& Encode(std::uint16_t const k) const override;
  std::uint16_t const& Decode(char16_t const code_point) const override;

 private:
  char16_t const* const blocks_begin_;
  char16_t const* const blocks_end_;
  std::array<char16_t, block_count * block_size> encoding_cache_;
  std::array<std::uint16_t, std::numeric_limits<char16_t>::max()>
      decoding_cache_;
};

template<std::size_t block_size, std::size_t block_count_plus_1>
constexpr CachingRepertoire<block_size, block_count_plus_1 - 1> MakeRepertoire(
    char16_t const (&blocks_begin)[block_count_plus_1]);

template<std::size_t block_size, std::size_t block_count>
template<std::size_t block_count_plus_1>
constexpr CachingRepertoire<block_size, block_count>::CachingRepertoire(
  char16_t const (&blocks_begin)[block_count_plus_1])
  : blocks_begin_(blocks_begin), blocks_end_(blocks_begin_ + block_count) {
  // Check null-termination.
  assert(blocks_begin_[block_count] == 0);
  // Check ordering and lack of overlap.
  for (char16_t const* block = blocks_begin_; block < blocks_end_; ++block) {
    //assert(*block < *(block + 1));
    //assert(*(block + 1) - *block >= block_size);
  }

  for (int block = 0; block < blocks_end_ - blocks_begin_; ++block) {
    char16_t const block_start_code_point = blocks_begin_[block];
    int const block_start_k = block_size * block;
    for (int offset = 0; offset < block_size; ++offset) {
      char16_t const code_point = block_start_code_point + offset;
      int const k = block_start_k + offset;
      encoding_cache_[k] = code_point;
      //decoding_cache_[code_point] = k;
    }
  }
}

template<std::size_t block_size, std::size_t block_count>
char16_t const& CachingRepertoire<block_size, block_count>::Encode(
    std::uint16_t const k) const {
  return encoding_cache_[k];
}

template<std::size_t block_size, std::size_t block_count>
std::uint16_t const& CachingRepertoire<block_size, block_count>::Decode(
    char16_t const code_point) const {
  return decoding_cache_[code_point];
}

template<std::size_t block_size, std::size_t block_count_plus_1>
constexpr CachingRepertoire<block_size, block_count_plus_1 - 1> MakeRepertoire(
    char16_t const (&blocks_begin)[block_count_plus_1]) {
  return CachingRepertoire<block_size, block_count_plus_1 - 1>(blocks_begin);
}

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

void Base32768Encode(Array<std::uint8_t const> input,
                     Array<std::uint8_t> output) {
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
    std::memcpy(
        output.data, &(repertoire->Encode(code_point)), sizeof(char16_t));

    // The following computation may cause |input.data| to overshoot the end if
    // using the special encoding at the end.  This is safe as soon as the loop
    // condition uses <.
    input_bit_index += bits_per_code_point;
    input.data += input_bit_index / bits_per_byte;
    input_bit_index %= bits_per_byte;
    output.data += sizeof(char16_t);
  }
}
UniqueArray<std::uint8_t> Base32768Encode(Array<std::uint8_t const> input,
                                          bool const null_terminated) {
  // TODO(phl): Add a function to compute the output size.
  base::UniqueArray<std::uint8_t> output((input.size << 1) +
                                         (null_terminated ? 1 : 0));
  if (output.size > 0) {
    base::Base32768Encode(input, output.get());
  }
  if (null_terminated) {
    output.data[output.size - 1] = 0;
  }
  return output;
}

void Base32768Decode(Array<std::uint8_t const> input,
                     Array<std::uint8_t> output) {
  CHECK_NOTNULL(input.data);
  CHECK_NOTNULL(output.data);
}

UniqueArray<std::uint8_t> Base32768Decode(Array<std::uint8_t const> input) {
  UniqueArray<std::uint8_t> output(input.size >> 1);
  if (output.size > 0) {
    Base32768Decode({ input.data, input.size & ~1 }, output.get());
  }
  return output;
}

}  // namespace internal_base32768
}  // namespace base
}  // namespace principia
