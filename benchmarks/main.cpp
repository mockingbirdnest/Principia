#include "benchmark/benchmark.h"
#include "absl/flags/parse.h"
#include "absl/log/check.h"
#include "absl/log/initialize.h"
#include "absl/log/log.h"

int __cdecl main(int argc, char* argv[]) {
  google::SetLogFilenameExtension(".log");
  absl::ParseCommandLine(argc, argv);
  absl::InitializeLog();
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
