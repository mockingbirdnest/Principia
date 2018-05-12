
// .\Release\x64\benchmarks.exe --benchmark_min_time=2 --benchmark_repetitions=10 --benchmark_filter=Base32768  // NOLINT(whitespace/line_length)

#include "base/hexadecimal.hpp"

#include <random>

#include "base/array.hpp"
#include "benchmark/benchmark.h"

namespace principia {
namespace base {

void BM_HexadecimalEncode(benchmark::State& state) {
  constexpr int min_input_size = 60'000;
  constexpr int max_input_size = 70'000;
  std::mt19937_64 random(42);
  std::uniform_int_distribution<std::uint64_t> size_distribution(
      min_input_size, max_input_size);
  std::uniform_int_distribution<int> bytes_distribution(0, 256);

  std::int64_t bytes_processed = 0;
  for (auto _ : state) {
    state.PauseTiming();
    UniqueArray<std::uint8_t> binary(size_distribution(random));
    for (int i = 0; i < binary.size; ++i) {
      binary.data[i] = bytes_distribution(random);
    }
    bytes_processed += binary.size;
    state.ResumeTiming();

    UniqueArray<char> const hexadecimal =
        HexadecimalEncode(binary.get(),
                          /*null_terminated=*/false);
    benchmark::DoNotOptimize(hexadecimal);
  }
  state.SetBytesProcessed(bytes_processed);
}

void BM_HexadecimalDecode(benchmark::State& state) {
  constexpr int preallocated_size = 1 << 20;
  constexpr int min_input_size = 120'000;
  constexpr int max_input_size = 140'000;

  std::mt19937_64 random(42);
  std::uniform_int_distribution<int> bytes_distribution(0, 256);

  // We need correct input data for the decoder.  Create it by enconding a large
  // chunk of data.
  UniqueArray<std::uint8_t> preallocated_binary(preallocated_size);
  for (int i = 0; i < preallocated_binary.size; ++i) {
    preallocated_binary.data[i] = bytes_distribution(random);
  }
  UniqueArray<char> const preallocated_hexadecimal =
      HexadecimalEncode(preallocated_binary.get(),
                      /*null_terminated=*/false);

  std::uniform_int_distribution<std::uint64_t> start_distribution(
      0, preallocated_hexadecimal.size - max_input_size);
  std::uniform_int_distribution<std::uint64_t> size_distribution(
      min_input_size, max_input_size);

  std::int64_t bytes_processed = 0;
  for (auto _ : state) {
    state.PauseTiming();
    auto const start = start_distribution(random);
    auto const size = size_distribution(random);
    Array<char> const hexadecimal(&preallocated_hexadecimal.data[start], size);
    bytes_processed += size * sizeof(char);
    state.ResumeTiming();

    UniqueArray<std::uint8_t> const binary = HexadecimalDecode(hexadecimal);
    benchmark::DoNotOptimize(binary);
  }
  state.SetBytesProcessed(bytes_processed);
}

BENCHMARK(BM_HexadecimalEncode);
BENCHMARK(BM_HexadecimalDecode);

}  // namespace base
}  // namespace principia
