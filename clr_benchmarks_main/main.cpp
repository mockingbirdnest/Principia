#include "clr_benchmarks/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "benchmark/benchmark.h"

int main(array<System::String^>^ args) {
  int argc = args->Length;
  char** argv = new char*[argc];
  for (int i = 0; i < argc; ++i) {
    array<System::Byte>^ encodedBytes =
        System::Text::Encoding::UTF8->GetBytes(args[i]);
    pin_ptr<System::Byte> pinnedBytes =
        &encodedBytes[encodedBytes->GetLowerBound(0)];
    argv[i] = new char[encodedBytes->Length + 1]; 
    std::memcpy(argv[i], reinterpret_cast<char*>(pinnedBytes),
                encodedBytes->Length);
    // null-terminate the native string
    argv[i][encodedBytes->Length] = '\0';
  }
  benchmark::Initialize(&argc, (char const**)(argv));
  principia::clr_benchmarks::SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator();

  benchmark::RunSpecifiedBenchmarks();
}
