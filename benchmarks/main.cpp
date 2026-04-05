#include "benchmark/benchmark.h"
#include "absl/log/check.h"
#include "absl/log/log.h"

int __cdecl main(int argc, char* argv[]) {
  google::SetLogFilenameExtension(".log");
  google::InitGoogleLogging(argv[0]);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
