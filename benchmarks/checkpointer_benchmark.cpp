// .\Release\x64\benchmarks.exe --benchmark_filter=Checkpointer --benchmark_repetitions=5  // NOLINT(whitespace/line_length)

#include "physics/checkpointer.hpp"

#include <memory>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {

using base::not_null;
using geometry::Instant;
using quantities::si::Second;
using serialization::Ephemeris;

// Stub writer for checkpointer. Does nothing.
void WriteToCheckpoint(not_null<typename Ephemeris::Checkpoint*> checkpoint) {}

// Stub reader for checkpointer. Does nothing.
absl::Status ReadFromCheckpoint(Ephemeris::Checkpoint const&) {
  return absl::OkStatus();
}

// Constructs a Checkpointer with the specified number of points. Checkpoints
// are spaced at an interval of one second.
std::unique_ptr<Checkpointer<Ephemeris>> ConstructCheckpointerWithSize(
    int size) {
  auto checkpointer = std::make_unique<Checkpointer<Ephemeris>>(
      &WriteToCheckpoint, &ReadFromCheckpoint);

  for (int i = 0; i < size; i++) {
    checkpointer->WriteToCheckpoint(Instant() + i * Second);
  }

  return checkpointer;
}

void BM_CheckpointerOldestCheckpoint(benchmark::State& state) {
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->oldest_checkpoint());
  }
}

void BM_CheckpointerNewestCheckpoint(benchmark::State& state) {
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->newest_checkpoint());
  }
}

void BM_CheckpointerCheckpointAtOrAfter(benchmark::State& state) {
  int const size = state.range(0);

  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(size);
  Instant const t = Instant() + (size / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->checkpoint_at_or_after(t));
  }
}
void BM_CheckpointerCheckpointAtOrBefore(benchmark::State& state) {
  int const size = state.range(0);

  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(size);
  Instant const t = Instant() + (size * 2.0 / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->checkpoint_at_or_before(t));
  }
}

void BM_CheckpointerAllCheckpoints(benchmark::State& state) {
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->all_checkpoints());
  }
}

void BM_CheckpointerAllCheckpointsAtOrBefore(benchmark::State& state) {
  int const size = state.range(0);
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(state.range(0));
  Instant const t = Instant() + (size / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->all_checkpoints_at_or_before(t));
  }
}
void BM_CheckpointerAllCheckpointsBetween(benchmark::State& state) {
  int const size = state.range(0);
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      ConstructCheckpointerWithSize(state.range(0));
  Instant const t1 = Instant() + (size / 3.0) * Second;
  Instant const t2 = Instant() + (size * 2.0 / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->all_checkpoints_between(t1, t2));
  }
}

BENCHMARK(BM_CheckpointerOldestCheckpoint)->Range(1, 512);
BENCHMARK(BM_CheckpointerNewestCheckpoint)->Range(1, 512);
BENCHMARK(BM_CheckpointerCheckpointAtOrAfter)->Range(1, 512);
BENCHMARK(BM_CheckpointerCheckpointAtOrBefore)->Range(1, 512);
BENCHMARK(BM_CheckpointerAllCheckpoints)->Range(1, 512);
BENCHMARK(BM_CheckpointerAllCheckpointsAtOrBefore)->Range(1, 512);
BENCHMARK(BM_CheckpointerAllCheckpointsBetween)->Range(1, 512);

}  // namespace physics
}  // namespace principia
