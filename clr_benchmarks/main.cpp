#include "benchmark/benchmark.h"

int main(int argc, const char* argv[]) {
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
