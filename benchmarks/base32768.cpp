
// .\Release\x64\benchmarks.exe --benchmark_min_time=2 --benchmark_repetitions=10 --benchmark_filter=Base32768  // NOLINT(whitespace/line_length)

#include "base/base32768.hpp"

#include <random>

#include "base/array.hpp"
#include "benchmark/benchmark.h"

namespace principia {
namespace base {

void BM_Base32768Encode(benchmark::State& state) {
  constexpr int min_input_size = 20'000;
  constexpr int max_input_size = 50'000;
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

    UniqueArray<char16_t> const base32768 =
        Base32768Encode(binary.get(),
                        /*null_terminated=*/false);
    benchmark::DoNotOptimize(base32768);
  }
  state.SetBytesProcessed(bytes_processed);
}

void BM_Base32768Decode(benchmark::State& state) {
  constexpr int preallocated_size = 1 << 20;
  constexpr int min_input_size = 10'000;
  constexpr int max_input_size = 25'000;

  std::mt19937_64 random(42);
  std::uniform_int_distribution<int> bytes_distribution(0, 256);

  // We need correct input data for the decoder.  Create it by enconding a large
  // chunk of data.
  UniqueArray<std::uint8_t> preallocated_binary(preallocated_size);
  for (int i = 0; i < preallocated_binary.size; ++i) {
    preallocated_binary.data[i] = bytes_distribution(random);
  }
  UniqueArray<char16_t> const preallocated_base32768 =
      Base32768Encode(preallocated_binary.get(),
                      /*null_terminated=*/false);

  std::uniform_int_distribution<std::uint64_t> start_distribution(
      0, preallocated_base32768.size - max_input_size);
  std::uniform_int_distribution<std::uint64_t> size_distribution(
      min_input_size, max_input_size);

  std::int64_t bytes_processed = 0;
  for (auto _ : state) {
    state.PauseTiming();
    auto const start = start_distribution(random);
    auto const size = size_distribution(random);
    Array<char16_t> const base32768(&preallocated_base32768.data[start], size);
    state.ResumeTiming();

    UniqueArray<std::uint8_t> const binary = Base32768Decode(base32768);
    bytes_processed += binary.size;
    benchmark::DoNotOptimize(binary);
  }
  state.SetBytesProcessed(bytes_processed);
}

BENCHMARK(BM_Base32768Encode);
BENCHMARK(BM_Base32768Decode);

}  // namespace base
}  // namespace principia
