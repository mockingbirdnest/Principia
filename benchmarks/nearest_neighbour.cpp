// .\Release\x64\benchmarks.exe --benchmark_filter=PCP --benchmark_repetitions=1  // NOLINT(whitespace/line_length)

#include "numerics/nearest_neighbour.hpp"

#include <random>
#include <vector>

#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace numerics {

using base::not_null;
using geometry::Frame;
using geometry::Vector;

using World = Frame<enum class WorldTag>;
using V = Vector<double, World>;

PrincipalComponentPartitioningTree<V> BuildTree(
    std::int64_t const points_in_tree,
    std::int64_t const max_values_per_cell,
    std::vector<V>& values) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate_distribution(-10, 10);

  std::vector<not_null<V const*>> pointers;
  values.reserve(points_in_tree);  // To avoid pointer invalidation below.
  values.clear();
  for (int i = 0; i < points_in_tree; ++i) {
    values.push_back(V({coordinate_distribution(random),
                        coordinate_distribution(random),
                        coordinate_distribution(random)}));
    pointers.push_back(&values.back());
  }
  return PrincipalComponentPartitioningTree<V>(pointers, max_values_per_cell);
}

void BM_PCPBuildTree(benchmark::State& state) {
  std::int64_t const points_in_tree = state.range(0);
  std::int64_t const max_values_per_cell = state.range(1);
  std::vector<V> values;

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        BuildTree(points_in_tree, max_values_per_cell, values));
  }
}

void BM_PCPFindNearestNeighbour(benchmark::State& state) {
  std::int64_t const points_in_tree = state.range(0);
  std::int64_t const max_values_per_cell = state.range(1);
  std::vector<V> values;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate_distribution(-10, 10);
  auto const tree = BuildTree(points_in_tree, max_values_per_cell, values);

  for (auto _ : state) {
    benchmark::DoNotOptimize(
        tree.FindNearestNeighbour(V({coordinate_distribution(random),
                                     coordinate_distribution(random),
                                     coordinate_distribution(random)})));
  }
}

BENCHMARK(BM_PCPBuildTree)
    ->Args({1'000, 1})
    ->Args({1'000, 4})
    ->Args({1'000, 16})
    ->Args({1'000, 64})
    ->Args({1'000, 256})
    ->Args({10'000, 1})
    ->Args({10'000, 4})
    ->Args({10'000, 16})
    ->Args({10'000, 64})
    ->Args({10'000, 256})
    ->Args({100'000, 1})
    ->Args({100'000, 4})
    ->Args({100'000, 16})
    ->Args({100'000, 64})
    ->Args({100'000, 256});
BENCHMARK(BM_PCPFindNearestNeighbour)
    ->Args({1'000, 1})
    ->Args({1'000, 4})
    ->Args({1'000, 16})
    ->Args({1'000, 64})
    ->Args({1'000, 256})
    ->Args({10'000, 1})
    ->Args({10'000, 4})
    ->Args({10'000, 16})
    ->Args({10'000, 64})
    ->Args({10'000, 256})
    ->Args({100'000, 1})
    ->Args({100'000, 4})
    ->Args({100'000, 16})
    ->Args({100'000, 64})
    ->Args({100'000, 256});

}  // namespace numerics
}  // namespace principia
