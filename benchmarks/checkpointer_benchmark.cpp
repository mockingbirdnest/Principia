// .\Release\x64\benchmarks.exe --benchmark_filter=Checkpointer --benchmark_repetitions=5  // NOLINT(whitespace/line_length)

#include "physics/checkpointer.hpp"

#include <memory>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "physics/checkpointer.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/integrators.pb.h"
#include "serialization/numerics.pb.h"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {

using base::not_null;
using geometry::Instant;
using quantities::si::Second;
using serialization::DoublePrecision;
using serialization::Ephemeris;
using serialization::IntegratorInstance;
using serialization::R3Element;
using serialization::SystemState;

// Stub writer for checkpointer. Adds some data to the checkpoint.
void WriteToCheckpoint(not_null<typename Ephemeris::Checkpoint*> checkpoint) {
  // Fill out the state with dummy values. We don't care what the values are (or
  // even if they're valid), just that they take time to copy.
  SystemState* state = checkpoint->mutable_instance()->mutable_current_state();
  state->mutable_time()->mutable_value()->set_double_(1.0);
  state->mutable_time()->mutable_error()->set_double_(1.0);
  for (int i = 0; i < 5; i++) {
    DoublePrecision* position = state->add_position();
    R3Element* raw = position->mutable_value()
                         ->mutable_point()
                         ->mutable_multivector()
                         ->mutable_vector();
    raw->mutable_x()->mutable_quantity()->set_magnitude(1.0);
    raw->mutable_y()->mutable_quantity()->set_magnitude(1.0);
    raw->mutable_y()->mutable_quantity()->set_magnitude(1.0);
    *state->add_velocity() = *position;
  }
}

// Stub reader for checkpointer. Does nothing.
absl::Status ReadFromCheckpoint(Ephemeris::Checkpoint const&) {
  return absl::OkStatus();
}

// Constructs a Checkpointer with the specified number of points. Checkpoints
// are spaced at an interval of one second.
std::unique_ptr<Checkpointer<Ephemeris>> NewCheckpointerWithSize(
    int const size) {
  auto checkpointer = std::make_unique<Checkpointer<Ephemeris>>(
      &WriteToCheckpoint, &ReadFromCheckpoint);

  for (int i = 0; i < size; i++) {
    checkpointer->WriteToCheckpoint(Instant() + i * Second);
  }

  return checkpointer;
}

void BM_CheckpointerOldestCheckpoint(benchmark::State& state) {
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->oldest_checkpoint());
  }
}

void BM_CheckpointerNewestCheckpoint(benchmark::State& state) {
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->newest_checkpoint());
  }
}

void BM_CheckpointerCheckpointAtOrAfter(benchmark::State& state) {
  int const size = state.range(0);
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(size);
  Instant const t = Instant() + (size / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->checkpoint_at_or_after(t));
  }
}
void BM_CheckpointerCheckpointAtOrBefore(benchmark::State& state) {
  int const size = state.range(0);
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(size);
  Instant const t = Instant() + (size * 2.0 / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->checkpoint_at_or_before(t));
  }
}

void BM_CheckpointerAllCheckpoints(benchmark::State& state) {
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->all_checkpoints());
  }
}

void BM_CheckpointerAllCheckpointsAtOrBefore(benchmark::State& state) {
  int const size = state.range(0);
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(state.range(0));
  Instant const t = Instant() + (size / 3.0) * Second;

  for (auto _ : state) {
    benchmark::DoNotOptimize(checkpointer->all_checkpoints_at_or_before(t));
  }
}
void BM_CheckpointerAllCheckpointsBetween(benchmark::State& state) {
  int const size = state.range(0);
  std::unique_ptr<Checkpointer<Ephemeris>> checkpointer =
      NewCheckpointerWithSize(state.range(0));
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
