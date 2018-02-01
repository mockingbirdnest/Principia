
#include "benchmark/benchmark.h"
#include "glog/logging.h"

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
