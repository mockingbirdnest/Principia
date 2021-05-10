
// .\Release\x64\benchmarks.exe --benchmark_filter=DiscreteTrajectory

#include "physics/discrete_trajectory.hpp"

#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "ksp_plugin/frames.hpp"

namespace principia {
namespace physics {

using base::make_not_null_unique;
using base::not_null;
using geometry::Displacement;
using geometry::Frame;
using geometry::Handedness;
using geometry::Inertial;
using geometry::Instant;
using geometry::Velocity;
using ksp_plugin::World;
using quantities::si::Metre;
using quantities::si::Second;

namespace {

// Creates a motionless trajectory with the given number of steps.
not_null<std::unique_ptr<DiscreteTrajectory<World>>> CreateMotionlessTrajectory(
    int const steps) {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  Instant t;
  for (int i = 0; i < steps; i++, t += 1 * Second) {
    trajectory->Append(t, {World::origin, World::unmoving});
  }
  return trajectory;
}

// Creates a circular trajectory with the given number of steps.
not_null<std::unique_ptr<DiscreteTrajectory<World>>> CreateCircularTrajectory(
    int const steps) {
  auto trajectory = make_not_null_unique<DiscreteTrajectory<World>>();
  Instant t;
  for (int i = 0; i < steps; i++, t += 1 * Second) {
    const double sine = std::sin(i);
    const double cosine = std::cos(i);
    trajectory->Append(
        t,
        {World::origin +
             Displacement<World>({cosine * Metre, 0 * Metre, sine * Metre}),
         Velocity<World>({-sine * Metre / Second,
                          0 * Metre / Second,
                          cosine * Metre / Second})});
  }
  return trajectory;
}

// Forks |parent| at a position |pos| of the way through.
// |parent| should be nonempty.
// |pos| should be in [0, 1].
not_null<DiscreteTrajectory<World>*> ForkAt(DiscreteTrajectory<World>& parent,
                                            double const pos) {
  CHECK(!parent.Empty());
  CHECK_GE(pos, 0);
  CHECK_LE(pos, 1);
  Instant const desired_fork_time =
      parent.t_min() + (parent.t_max() - parent.t_min()) * pos;
  auto const fork_it = parent.LowerBound(desired_fork_time);
  return parent.NewForkWithCopy(fork_it->time);
}

}  // namespace

void BM_DiscreteTrajectoryFront(benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(4);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    fork->front();
  }
}

void BM_DiscreteTrajectoryBack(benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(4);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    fork->back();
  }
}

void BM_DiscreteTrajectoryBegin(benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(4);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    fork->begin();
  }
}

void BM_DiscreteTrajectoryEnd(benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(4);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    fork->end();
  }
}

void BM_DiscreteTrajectoryTMin(benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(4);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    fork->t_min();
  }
}

void BM_DiscreteTrajectoryTMax(benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(4);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    fork->t_max();
  }
}

void BM_DiscreteTrajectoryCreateDestroy(benchmark::State& state) {
  int const steps = state.range(0);
  for (auto _ : state) {
    CreateMotionlessTrajectory(steps);
  }
}

void BM_DiscreteTrajectoryIterate(benchmark::State& state) {
  int const steps = state.range(0);
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(steps);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    for (auto it = fork->begin(); it != fork->end(); ++it) {
    }
  }
}

void BM_DiscreteTrajectoryReverseIterate(benchmark::State& state) {
  int const steps = state.range(0);
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(steps);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    for (auto it = fork->end(); it != fork->begin(); --it) {
    }
  }
}

void BM_DiscreteTrajectoryFind(benchmark::State& state) {
  int const steps = state.range(0);
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(steps);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    // These times are in different segments of the fork.
    fork->Find(Instant() + 1.0 / 3.0 * steps * Second);
    fork->Find(Instant() + 2.0 / 3.0 * steps * Second);
    fork->Find(Instant() + 5.0 / 6.0 * steps * Second);
  }
}

void BM_DiscreteTrajectoryLowerBound(benchmark::State& state) {
  int const steps = state.range(0);
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateMotionlessTrajectory(steps);
  not_null<DiscreteTrajectory<World>*> const fork =
      ForkAt(*ForkAt(*trajectory, 0.5), 0.75);

  for (auto _ : state) {
    // These times are in different segments of the fork.
    fork->LowerBound(Instant() + 1.0 / 3.0 * steps * Second);
    fork->LowerBound(Instant() + 2.0 / 3.0 * steps * Second);
    fork->LowerBound(Instant() + 5.0 / 6.0 * steps * Second);
  }
}

void BM_DiscreteTrajectoryEvaluateDegreesOfFreedomExact(
    benchmark::State& state) {
  int const steps = state.range(0);
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateCircularTrajectory(4);

  Instant const t = Instant() + 1 * Second;
  for (auto _ : state) {
    trajectory->EvaluateDegreesOfFreedom(t);
  }
}

void BM_DiscreteTrajectoryEvaluateDegreesOfFreedomInterpolated(
    benchmark::State& state) {
  not_null<std::unique_ptr<DiscreteTrajectory<World>>> const trajectory =
      CreateCircularTrajectory(4);

  Instant const t = Instant() + 1.5 * Second;
  for (auto _ : state) {
    trajectory->EvaluateDegreesOfFreedom(t);
  }
}

BENCHMARK(BM_DiscreteTrajectoryFront);
BENCHMARK(BM_DiscreteTrajectoryBack);
BENCHMARK(BM_DiscreteTrajectoryBegin);
BENCHMARK(BM_DiscreteTrajectoryEnd);
BENCHMARK(BM_DiscreteTrajectoryTMin);
BENCHMARK(BM_DiscreteTrajectoryTMax);
BENCHMARK(BM_DiscreteTrajectoryCreateDestroy)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryIterate)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryReverseIterate)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryFind)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryLowerBound)->Range(8, 1024);
BENCHMARK(BM_DiscreteTrajectoryEvaluateDegreesOfFreedomExact);
BENCHMARK(BM_DiscreteTrajectoryEvaluateDegreesOfFreedomInterpolated);

}  // namespace physics
}  // namespace principia
