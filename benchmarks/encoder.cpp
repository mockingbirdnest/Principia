
// .\Release\x64\benchmarks.exe --benchmark_min_time=2 --benchmark_repetitions=10 --benchmark_filter=Base32768  // NOLINT(whitespace/line_length)

#include "base/encoder.hpp"

#include <random>

#include "base/array.hpp"
#include "base/base32768.hpp"
#include "base/hexadecimal.hpp"
#include "benchmark/benchmark.h"

// Clang doesn't have a correct |std::array| yet, and we don't actually use this
// code, so let's get rid of the entire body.
#if PRINCIPIA_COMPILER_MSVC

namespace principia {
namespace base {

template<typename Encoder>
void BM_Encode(benchmark::State& state) {
  constexpr int preallocated_size = 1 << 20;
  constexpr int min_input_size = 20'000;
  constexpr int max_input_size = 50'000;

  Encoder encoder;
  std::mt19937_64 random(42);
  std::uniform_int_distribution<int> bytes_distribution(0, 256);

  UniqueArray<std::uint8_t> preallocated_binary(preallocated_size);
  for (int i = 0; i < preallocated_binary.size; ++i) {
    preallocated_binary.data[i] = bytes_distribution(random);
  }

  std::uniform_int_distribution<std::uint64_t> start_distribution(
      0, preallocated_binary.size - max_input_size);
  std::uniform_int_distribution<std::uint64_t> size_distribution(
      min_input_size, max_input_size);

  std::int64_t bytes_processed = 0;
  for (auto _ : state) {
    state.PauseTiming();
    auto const start = start_distribution(random);
    auto const size = size_distribution(random);
    Array<std::uint8_t> binary(&preallocated_binary.data[start], size);
    bytes_processed += binary.size;
    state.ResumeTiming();

    UniqueArray<typename Encoder::Char> const encoded = encoder.Encode(binary);
    benchmark::DoNotOptimize(encoded);
  }
  state.SetBytesProcessed(bytes_processed);
}

template<typename Encoder>
void BM_Decode(benchmark::State& state) {
  constexpr int preallocated_size = 1 << 20;
  constexpr int min_input_size = 10'000;
  constexpr int max_input_size = 25'000;

  Encoder encoder;
  std::mt19937_64 random(42);
  std::uniform_int_distribution<int> bytes_distribution(0, 256);

  // We need correct input data for the decoder.  Create it by enconding a large
  // chunk of data.
  UniqueArray<std::uint8_t> preallocated_binary(preallocated_size);
  for (int i = 0; i < preallocated_binary.size; ++i) {
    preallocated_binary.data[i] = bytes_distribution(random);
  }
  UniqueArray<typename Encoder::Char> const preallocated_encoded =
      encoder.Encode(preallocated_binary.get());

  std::uniform_int_distribution<std::uint64_t> start_distribution(
      0, preallocated_encoded.size - max_input_size);
  std::uniform_int_distribution<std::uint64_t> size_distribution(
      min_input_size, max_input_size);

  std::int64_t bytes_processed = 0;
  for (auto _ : state) {
    state.PauseTiming();
    auto const start = start_distribution(random);
    auto const size = size_distribution(random);
    Array<typename Encoder::Char> const encoded(
        &preallocated_encoded.data[start], size);
    state.ResumeTiming();

    UniqueArray<std::uint8_t> const binary = encoder.Decode(encoded);
    bytes_processed += binary.size;
    benchmark::DoNotOptimize(binary);
  }
  state.SetBytesProcessed(bytes_processed);
}

using Encoder16 = HexadecimalEncoder</*null_terminated=*/false>;
using Encoder32768 = Base32768Encoder</*null_terminated=*/false>;

BENCHMARK_TEMPLATE(BM_Encode, Encoder16);
BENCHMARK_TEMPLATE(BM_Decode, Encoder16);
#if !PRINCIPIA_COMPILER_MSVC || \
    !(_MSC_FULL_VER == 191526608 || \
      _MSC_FULL_VER == 191526731 || \
      _MSC_FULL_VER == 191627024)
BENCHMARK_TEMPLATE(BM_Encode, Encoder32768);
BENCHMARK_TEMPLATE(BM_Decode, Encoder32768);
#endif

}  // namespace base
}  // namespace principia

#endif
