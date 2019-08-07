
#include "benchmark/benchmark.h"
#include "glog/logging.h"

int main(int argc, char* argv[]) {
  google::SetLogFilenameExtension(".log");
  google::InitGoogleLogging(argv[0]);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
