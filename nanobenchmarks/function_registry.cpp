#include "nanobenchmarks/function_registry.hpp"

#include <functional>
#include <map>
#include <string>

#include "glog/logging.h"

namespace principia {
namespace nanobenchmarks {
namespace _function_registry {
namespace internal {

std::string const& Nanobenchmark::name() const {
  return name_;
}

void Nanobenchmark::SetName(std::string_view const name) {
  name_ = name;
}

LatencyDistributionTable Fixture::Run() {
  std::size_t const sample_count = absl::GetFlag(FLAGS_samples);
  std::size_t const loop_iterations = absl::GetFlag(FLAGS_loop_iterations);
  static std::vector<double>& samples = *new std::vector<double>(
      sample_count, std::numeric_limits<double>::quiet_NaN());
  int registers[4]{};
  int leaf = 0;
  for (int j = 0; j < sample_count; ++j) {
    double const input = absl::GetFlag(FLAGS_input);
    double x = input;
    // The CPUID barriers prevent out-of-order execution; see [Pao10].
#if PRINCIPIA_COMPILER_MSVC
    __cpuid(registers, leaf);
#else
    asm volatile("cpuid" ::: "eax", "ebx", "ecx", "edx", "memory");
#endif
    auto const tsc_start = __rdtsc();
    for (int i = 0; i < loop_iterations; ++i) {
      x = NanobenchmarkCase(x);
      x += input - x;
    }
    unsigned int tsc_aux;
    // The use of RDTSCP rather than RDTSC here follows [Pao10].  See the
    // Intel® 64 and IA-32 Architectures Software Developer’s Manual: The
    // RDTSCP  instruction is not a serializing instruction, but it does wait
    // until all previous instructions have executed and all previous loads
    // are globally visible.  But it does not wait for previous stores to be
    // globally visible, and subsequent instructions may begin execution
    // before the read operation is performed.
    auto const tsc_stop = __rdtscp(&tsc_aux);
#if PRINCIPIA_COMPILER_MSVC
    __cpuid(registers, leaf);
#else
    asm volatile("cpuid" ::: "eax", "ebx", "ecx", "edx", "memory");
#endif
    double const δtsc = tsc_stop - tsc_start;
    samples[j] = δtsc / loop_iterations;
  }
  // if (logger != nullptr) {
  //   logger->Append(
  //       "samples",
  //       std::tuple{FunctionRegistry::names_by_function().at(f), samples});
  // }
  std::ranges::sort(samples);
  LatencyDistributionTable result(absl::GetFlag(FLAGS_quantiles));
  result.SetSamples(samples);
  return result;
}

bool FunctionRegistry::Register(std::string_view name,
                                BenchmarkedFunction function) {
  CHECK(singleton().names_by_function_.emplace(function, name).second)
      << " Registering function " << function << " as " << name << ": "
      << "function already registered as "
      << singleton().names_by_function_[function];
  CHECK(singleton().functions_by_name_.emplace(name, function).second)
      << " Registering function " << function << " as " << name << ": "
      << " name already taken by "
      << singleton().functions_by_name_.find(name)->second;
  return true;
}

FunctionRegistry& FunctionRegistry::singleton() {
  static auto& singleton = *new FunctionRegistry;
  return singleton;
}

std::map<std::string, BenchmarkedFunction, std::less<>> const&
FunctionRegistry::functions_by_name() {
  return singleton().functions_by_name_;
}

std::map<BenchmarkedFunction, std::string> const&
FunctionRegistry::names_by_function() {
  return singleton().names_by_function_;
}

Nanobenchmark* NanobenchmarkRegistry::Register(
    Nanobenchmark* const nanobenchmark) {
  singleton().nanobenchmarks_by_name_.emplace(nanobenchmark->name(),
                                              nanobenchmark);
  return nanobenchmark;
}

std::vector<Nanobenchmark*> const&
NanobenchmarkRegistry::NanobenchmarksMatching(std::regex const& filter) {
  std::vector<Nanobenchmark*> matching;
  for (auto const& [name, nanobenchmark] :
       singleton().nanobenchmarks_by_name_) {
    if (std::regex_match(name, filter)) {
      matching.push_back(nanobenchmark);
    }
  }
  return matching;
}

}  // namespace internal
}  // namespace _function_registry
}  // namespace nanobenchmarks
}  // namespace principia
