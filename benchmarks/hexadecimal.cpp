
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include "base/hexadecimal.hpp"
#include "base/not_null.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

using principia::base::not_null;
using principia::base::HexadecimalDecode;
using principia::base::HexadecimalEncode;

namespace principia {
namespace benchmarks {

static char const kPi1000HexadecimalDigits[] =
    "3243F6A8885A308D313198A2E03707344A4093822299F31D0082EFA98EC4E6C89452821E63"
    "8D01377BE5466CF34E90C6CC0AC29B7C97C50DD3F84D5B5B54709179216D5D98979FB1BD13"
    "10BA698DFB5AC2FFD72DBD01ADFB7B8E1AFED6A267E96BA7C9045F12C7F9924A19947B3916"
    "CF70801F2E2858EFC16636920D871574E69A458FEA3F4933D7E0D95748F728EB658718BCD5"
    "882154AEE7B54A41DC25A59B59C30D5392AF26013C5D1B023286085F0CA417918B8DB38EF8"
    "E79DCB0603A180E6C9E0E8BB01E8A3ED71577C1BD314B2778AF2FDA55605C60E65525F3AA5"
    "5AB945748986263E8144055CA396A2AAB10B6B4CC5C341141E8CEA15486AF7C72E993B3EE1"
    "411636FBC2A2BA9C55D741831F6CE5C3E169B87931EAFD6BA336C24CF5C7A3253812895867"
    "73B8F48986B4BB9AFC4BFE81B6628219361D809CCFB21A991487CAC605DEC8032EF845D5DE"
    "98575B1DC262302EB651B8823893E81D396ACC50F6D6FF383F442392E0B4482A484200469C"
    "8F04A9E1F9B5E21C66842F6E96C9A670C9C61ABD388F06A51A0D2D8542F68960FA728AB513"
    "3A36EEF0B6C137A3BE4BA3BF0507EFB2A98A1F1651D39AF017666CA593E82430E888CEE861"
    "9456F9FB47D84A5C33B8B5EBEE06F75D885C12073401A449F56C16AA64ED3AA62363F77061"
    "BFEDF72429B023D37D0D724D00A1248DB0FEAD";

static char const kPi500Bytes[] =
    "\x32\x43\xF6\xA8\x88\x5A\x30\x8D\x31\x31\x98\xA2\xE0\x37\x07\x34\x4A\x40"
    "\x93\x82\x22\x99\xF3\x1D\x00\x82\xEF\xA9\x8E\xC4\xE6\xC8\x94\x52\x82\x1E"
    "\x63\x8D\x01\x37\x7B\xE5\x46\x6C\xF3\x4E\x90\xC6\xCC\x0A\xC2\x9B\x7C\x97"
    "\xC5\x0D\xD3\xF8\x4D\x5B\x5B\x54\x70\x91\x79\x21\x6D\x5D\x98\x97\x9F\xB1"
    "\xBD\x13\x10\xBA\x69\x8D\xFB\x5A\xC2\xFF\xD7\x2D\xBD\x01\xAD\xFB\x7B\x8E"
    "\x1A\xFE\xD6\xA2\x67\xE9\x6B\xA7\xC9\x04\x5F\x12\xC7\xF9\x92\x4A\x19\x94"
    "\x7B\x39\x16\xCF\x70\x80\x1F\x2E\x28\x58\xEF\xC1\x66\x36\x92\x0D\x87\x15"
    "\x74\xE6\x9A\x45\x8F\xEA\x3F\x49\x33\xD7\xE0\xD9\x57\x48\xF7\x28\xEB\x65"
    "\x87\x18\xBC\xD5\x88\x21\x54\xAE\xE7\xB5\x4A\x41\xDC\x25\xA5\x9B\x59\xC3"
    "\x0D\x53\x92\xAF\x26\x01\x3C\x5D\x1B\x02\x32\x86\x08\x5F\x0C\xA4\x17\x91"
    "\x8B\x8D\xB3\x8E\xF8\xE7\x9D\xCB\x06\x03\xA1\x80\xE6\xC9\xE0\xE8\xBB\x01"
    "\xE8\xA3\xED\x71\x57\x7C\x1B\xD3\x14\xB2\x77\x8A\xF2\xFD\xA5\x56\x05\xC6"
    "\x0E\x65\x52\x5F\x3A\xA5\x5A\xB9\x45\x74\x89\x86\x26\x3E\x81\x44\x05\x5C"
    "\xA3\x96\xA2\xAA\xB1\x0B\x6B\x4C\xC5\xC3\x41\x14\x1E\x8C\xEA\x15\x48\x6A"
    "\xF7\xC7\x2E\x99\x3B\x3E\xE1\x41\x16\x36\xFB\xC2\xA2\xBA\x9C\x55\xD7\x41"
    "\x83\x1F\x6C\xE5\xC3\xE1\x69\xB8\x79\x31\xEA\xFD\x6B\xA3\x36\xC2\x4C\xF5"
    "\xC7\xA3\x25\x38\x12\x89\x58\x67\x73\xB8\xF4\x89\x86\xB4\xBB\x9A\xFC\x4B"
    "\xFE\x81\xB6\x62\x82\x19\x36\x1D\x80\x9C\xCF\xB2\x1A\x99\x14\x87\xCA\xC6"
    "\x05\xDE\xC8\x03\x2E\xF8\x45\xD5\xDE\x98\x57\x5B\x1D\xC2\x62\x30\x2E\xB6"
    "\x51\xB8\x82\x38\x93\xE8\x1D\x39\x6A\xCC\x50\xF6\xD6\xFF\x38\x3F\x44\x23"
    "\x92\xE0\xB4\x48\x2A\x48\x42\x00\x46\x9C\x8F\x04\xA9\xE1\xF9\xB5\xE2\x1C"
    "\x66\x84\x2F\x6E\x96\xC9\xA6\x70\xC9\xC6\x1A\xBD\x38\x8F\x06\xA5\x1A\x0D"
    "\x2D\x85\x42\xF6\x89\x60\xFA\x72\x8A\xB5\x13\x3A\x36\xEE\xF0\xB6\xC1\x37"
    "\xA3\xBE\x4B\xA3\xBF\x05\x07\xEF\xB2\xA9\x8A\x1F\x16\x51\xD3\x9A\xF0\x17"
    "\x66\x6C\xA5\x93\xE8\x24\x30\xE8\x88\xCE\xE8\x61\x94\x56\xF9\xFB\x47\xD8"
    "\x4A\x5C\x33\xB8\xB5\xEB\xEE\x06\xF7\x5D\x88\x5C\x12\x07\x34\x01\xA4\x49"
    "\xF5\x6C\x16\xAA\x64\xED\x3A\xA6\x23\x63\xF7\x70\x61\xBF\xED\xF7\x24\x29"
    "\xB0\x23\xD3\x7D\x0D\x72\x4D\x00\xA1\x24\x8D\xB0\xFE\xAD";

int const kCopiesOfPi = 10000;

template<typename Container>
std::enable_if_t<
    std::is_convertible<typename Container::value_type, uint8_t>::value,
    void>
HexadecimalEncodeNoMemcpy(Container const& input, not_null<Container*> output) {
  Container const& bytes = input;
  Container digits;
  digits.resize(bytes.size() * 2);
  // The following was undefined behaviour pre-C++11, but it is now well-defined
  // even when |digits.size() == 0|.  We do not use |digits.data()| because this
  // only works for |std::vector| (it is read-only in a |std::basic_string|).
  char* digit = &digits[0];
  for (uint8_t const byte : bytes) {
    auto blurp = byte << 1;
    *digit = base::kByteToHexadecimalDigits[blurp];
    *++digit = base::kByteToHexadecimalDigits[++blurp];
    ++digit;
  }
  *output = std::move(digits);
}

void HexEncodeNoMemcpy(not_null<benchmark::State*> const state,
                      not_null<bool*> const correct,
                      std::string const& input_bytes,
                      std::string const& expected_digits,
                      std::string* const digits) {
  HexadecimalEncodeNoMemcpy<std::string>(input_bytes, digits);
  state->PauseTiming();
  *correct = *digits == expected_digits;
  state->ResumeTiming();
}

void BM_EncodePiNoMemcpy(benchmark::State& state) {  // NOLINT(runtime/references)
  bool correct;
  state.PauseTiming();
  std::string const pi_bytes(kPi500Bytes, 500);
  std::string input_bytes;
  std::string expected_digits;
  input_bytes.reserve(500 * kCopiesOfPi);
  expected_digits.reserve(1000 * kCopiesOfPi);
  for (int i = kCopiesOfPi; i > 0; --i) {
    input_bytes += pi_bytes;
    expected_digits += kPi1000HexadecimalDigits;
  }
  state.ResumeTiming();
  while (state.KeepRunning()) {
    state.PauseTiming();
    std::string buffer;
    state.ResumeTiming();
    HexEncodeNoMemcpy(&state, &correct, input_bytes, expected_digits, &buffer);
  }
  std::stringstream ss;
  ss << correct;
  state.SetLabel(ss.str());
}

BENCHMARK(BM_EncodePiNoMemcpy);

void HexEncode(not_null<benchmark::State*> const state,
              not_null<bool*> const correct,
              std::string const& input_bytes,
              std::string const& expected_digits,
              std::string* const digits) {
  HexadecimalEncode<std::string>(input_bytes, digits);
  state->PauseTiming();
  *correct = *digits == expected_digits;
  state->ResumeTiming();
}

void BM_EncodePi(benchmark::State& state) {  // NOLINT(runtime/references)
  bool correct;
  state.PauseTiming();
  std::string const pi_bytes(kPi500Bytes, 500);
  std::string input_bytes;
  std::string expected_digits;
  input_bytes.reserve(500 * kCopiesOfPi);
  expected_digits.reserve(1000 * kCopiesOfPi);
  for (int i = kCopiesOfPi; i > 0; --i) {
    input_bytes += pi_bytes;
    expected_digits += kPi1000HexadecimalDigits;
  }
  state.ResumeTiming();
  while (state.KeepRunning()) {
    state.PauseTiming();
    std::string buffer;
    state.ResumeTiming();
    HexEncode(&state, &correct, input_bytes, expected_digits, &buffer);
  }
  std::stringstream ss;
  ss << correct;
  state.SetLabel(ss.str());
}

BENCHMARK(BM_EncodePi);

void BM_EncodePiNoMemcpyPreallocated(
    benchmark::State& state) {  // NOLINT(runtime/references)
  bool correct;
  state.PauseTiming();
  std::string const pi_bytes(kPi500Bytes, 500);
  std::string input_bytes;
  std::string expected_digits;
  input_bytes.reserve(500 * kCopiesOfPi);
  expected_digits.reserve(1000 * kCopiesOfPi);
  for (int i = kCopiesOfPi; i > 0; --i) {
    input_bytes += pi_bytes;
    expected_digits += kPi1000HexadecimalDigits;
  }
  std::string buffer;
  buffer.resize(expected_digits.size());
  state.ResumeTiming();
  while (state.KeepRunning()) {
    HexEncodeNoMemcpy(&state, &correct, input_bytes, expected_digits, &buffer);
  }
  std::stringstream ss;
  ss << correct;
  state.SetLabel(ss.str());
}

BENCHMARK(BM_EncodePiNoMemcpyPreallocated);

void BM_EncodePiPreallocated(
    benchmark::State& state) {  // NOLINT(runtime/references)
  bool correct;
  state.PauseTiming();
  std::string const pi_bytes(kPi500Bytes, 500);
  std::string input_bytes;
  std::string expected_digits;
  input_bytes.reserve(500 * kCopiesOfPi);
  expected_digits.reserve(1000 * kCopiesOfPi);
  for (int i = kCopiesOfPi; i > 0; --i) {
    input_bytes += pi_bytes;
    expected_digits += kPi1000HexadecimalDigits;
  }
  std::string buffer;
  buffer.resize(expected_digits.size());
  state.ResumeTiming();
  while (state.KeepRunning()) {
    HexEncode(&state, &correct, input_bytes, expected_digits, &buffer);
  }
  std::stringstream ss;
  ss << correct;
  state.SetLabel(ss.str());
}

BENCHMARK(BM_EncodePiPreallocated);

}  // namespace benchmarks
}  // namespace principia
