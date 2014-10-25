#include "benchmark/benchmark.h"

int main(int argc, char const* argv[]) {
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
