// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_filter=ThreadPool  // NOLINT(whitespace/line_length)

#include <cmath>
#include <cstdint>
#include <future>
#include <random>
#include <vector>

#include "absl/synchronization/mutex.h"
#include "base/thread_pool.hpp"
#include "benchmark/benchmark.h"

namespace principia {
namespace base {

using namespace principia::base::_thread_pool;

absl::Mutex lock;
std::mt19937_64 random(42);
std::uniform_int_distribution<int> distribution(0, 1e5);

double ComsumeCpuNoLock(std::int64_t const n) {
  double result;
  {
    absl::ReaderMutexLock l(&lock);
    result = distribution(random);
  }
  for (int i = 0; i < n; ++i) {
    result += std::sqrt(i);
  }
  return result;
}

double ComsumeCpuSharedLock(std::int64_t const n) {
  absl::ReaderMutexLock l(&lock);
  double result = distribution(random);
  for (int i = 0; i < n; ++i) {
    result += std::sqrt(i);
  }
  return result;
}

double ComsumeCpuExclusiveLock(std::int64_t const n) {
  absl::MutexLock l(&lock);
  double result = distribution(random);
  for (int i = 0; i < n; ++i) {
    result += std::sqrt(i);
  }
  return result;
}

void BM_ThreadPoolNoLock(benchmark::State& state) {
  ThreadPool<void> pool(/*pool_size=*/state.range(0));
  std::vector<std::int64_t> results;
  for (auto _ : state) {
    std::vector<std::future<void>> futures;
    for (int i = 0; i < 1000; ++i) {
      futures.push_back(pool.Add([]() {
        double const result = ComsumeCpuNoLock(1e5);
        benchmark::DoNotOptimize(result);
      }));
    }
    for (auto const& future : futures) {
      future.wait();
    }
  }
}

void BM_ThreadPoolSharedLock(benchmark::State& state) {
  ThreadPool<void> pool(/*pool_size=*/state.range(0));
  std::vector<std::int64_t> results;
  for (auto _ : state) {
    std::vector<std::future<void>> futures;
    for (int i = 0; i < 1000; ++i) {
      futures.push_back(pool.Add([]() {
        double const result = ComsumeCpuSharedLock(1e5);
        benchmark::DoNotOptimize(result);
      }));
    }
    for (auto const& future : futures) {
      future.wait();
    }
  }
}

void BM_ThreadPoolExclusiveLock(benchmark::State& state) {
  ThreadPool<void> pool(/*pool_size=*/state.range(0));
  std::vector<std::int64_t> results;
  for (auto _ : state) {
    std::vector<std::future<void>> futures;
    for (int i = 0; i < 1000; ++i) {
      futures.push_back(pool.Add([]() {
        double const result = ComsumeCpuExclusiveLock(1e5);
        benchmark::DoNotOptimize(result);
      }));
    }
    for (auto const& future : futures) {
      future.wait();
    }
  }
}

BENCHMARK(BM_ThreadPoolNoLock)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Arg(6)
    ->Arg(7)
    ->Arg(8)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_ThreadPoolSharedLock)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Arg(6)
    ->Arg(7)
    ->Arg(8)
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_ThreadPoolExclusiveLock)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Arg(6)
    ->Arg(7)
    ->Arg(8)
    ->Unit(benchmark::kMillisecond);

}  // namespace base
}  // namespace principia
