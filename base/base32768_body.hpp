
#pragma once

#include "base/base32768.hpp"

#include <array>
#include <cstdint>
#include <cstring>
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

class Repertoire {
 public:
  virtual ~Repertoire() = default;

  virtual char16_t const& Encode(std::uint16_t const k) = 0;
  virtual std::uint16_t const& Decode(char16_t const code_point) = 0;
};

template<std::size_t block_count, std::size_t block_size>
class CachingRepertoire : public Repertoire {
 public:
  char16_t const& Encode(std::uint16_t const k) override;
  std::uint16_t const& Decode(char16_t const code_point) override;

 private:
  template<std::size_t block_count_plus_1>
  constexpr CachingRepertoire(
    char16_t const (&begin_block)[block_count_plus_1]);

  char16_t const* begin_block_;
  std::array<char16_t, block_count * block_size> encoding_cache_;
  std::array<std::uint16_t, sizeof(char16_t)> decoding_cache_;
};

template<std::size_t block_count, std::size_t block_size>
template<std::size_t block_count_plus_1>
constexpr CachingRepertoire<block_count, block_size>::CachingRepertoire(
    char16_t const (&begin_block)[block_count_plus_1])
    : begin_block_(begin_block) {
  // Check null-termination.
  assert(begin_block_[block_count] == 0);
  // Check ordering and lack of overlap.
  for (char16_t const* block = begin_block_;
       block < begin_block_ + block_count_ - 1;
       ++block) {
    assert(*block < *(block + 1));
    assert(*(block + 1) - *block >= block_size);
  }
}

constexpr std::array<char16_t const*, 2> repertoire = {
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
    u"ꊠꋀꋠꌀꌠꍀꍠꎀꎠꏀꏠꐀꐠꑀꑠ꒠ꔀꔠꕀꕠꖀꖠꗀꗠꙀꚠꛀ꜀꜠ꝀꞀꡀ",  // length = 1 << 10
    u"ƀɀɠʀ",        // length = 1 << 2
};

constexpr int block_size = 1 << 5;

char16_t const& Encode(Repertoire const r, int const k) {
  static std::array<std::map<int, char16_t>*, 2> encode = []() {
    std::array<std::map<int, char16_t>*, 2> result;
    for (int r = Begin; r < End; ++r) {
      result[r] = new std::map<int, char16_t>;
      for (int block = 0;
        block < std::char_traits<char16_t>::length(repertoire[r]);
        ++block) {
        char16_t const block_start_code_point = repertoire[r][block];
        int const block_start_k = block_size * block;
        for (int offset = 0; offset < block_size; ++offset) {
          char16_t const code_point = block_start_code_point + offset;
          int const k = block_start_k + offset;
          (*result[r])[k] = code_point;
        }
      }
    }
    return result;
  }();
  return (*encode[r])[k];
}

int const& Decode(Repertoire const r, char16_t const code_point) {
  static std::array<std::map<char16_t, int>*, 2> decode = []() {
    std::array<std::map<char16_t, int>*, 2> result;
    for (int r = Begin; r < End; ++r) {
      result[r] = new std::map<char16_t, int>;
      for (int block = 0;
        block < std::char_traits<char16_t>::length(repertoire[r]);
        ++block) {
        char16_t const block_start_code_point = repertoire[r][block];
        int const block_start_k = block_size * block;
        for (int offset = 0; offset < block_size; ++offset) {
          char16_t const code_point = block_start_code_point + offset;
          int const k = block_start_k + offset;
          (*result[r])[code_point] = k;
        }
      }
    }
    return result;
  }();
  return (*decode[r])[code_point];
}

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
    Repertoire repertoire = FifteenBits;

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
        repertoire = SevenBits;
      }
    }
    std::int32_t code_point = (data & mask) >> shift;
    CHECK_LE(0, code_point);
    CHECK_LT(code_point, 1 << bits_per_code_point);
    std::memcpy(output.data, &Encode(repertoire, code_point), sizeof(char16_t));

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
