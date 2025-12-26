#pragma once

#include "nanobenchmarks/nanobenchmark.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "glog/logging.h"
#include "nanobenchmarks/dependencies.hpp"

#if PRINCIPIA_COMPILER_MSVC
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

namespace principia {
namespace nanobenchmarks {
namespace _nanobenchmark {
namespace internal {

using namespace principia::nanobenchmarks::_dependencies;

template<typename Value_, typename Argument_>
Nanobenchmark<Value_, Argument_>::BenchmarkedFunction
Nanobenchmark<Value_, Argument_>::function() const {
  return function_;
}

template<typename Value_, typename Argument_>
std::string const& Nanobenchmark<Value_, Argument_>::name() const {
  return name_;
}

template<typename Value_, typename Argument_>
void Nanobenchmark<Value_, Argument_>::SetFunction(
    BenchmarkedFunction function) {
  CHECK(function != nullptr);
  function_ = function;
}

template<typename Value_, typename Argument_>
void Nanobenchmark<Value_, Argument_>::SetName(std::string_view const name) {
  name_ = name;
}

template<typename Value_, typename Argument_>
__declspec(noinline) LatencyDistributionTable
Nanobenchmark<Value_, Argument_>::Run(Logger* const logger) const {
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
      using D = Dependencies<Value, Argument>;
      x = D::ConsumeValue(NanobenchmarkCase(D::ProduceArgument(x)));
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
  if (logger != nullptr) {
    logger->Append("samples", std::tuple{name(), samples});
  }
  std::ranges::sort(samples);
  LatencyDistributionTable result;
  result.SetSamples(samples);
  return result;
}

template<typename Value_, typename Argument_>
Nanobenchmark<Value_, Argument_> const*
NanobenchmarkRegistry<Value_, Argument_>::Register(
    Nanobenchmark<Value_, Argument_> const* const nanobenchmark) {
  singleton().nanobenchmarks_by_name_.emplace(nanobenchmark->name(),
                                              nanobenchmark);
  if (nanobenchmark->function() != nullptr) {
    singleton().nanobenchmarks_by_function_.emplace(nanobenchmark->function(),
                                                    nanobenchmark);
  }
  return nanobenchmark;
}

template<typename Value_, typename Argument_>
std::vector<Nanobenchmark<Value_, Argument_> const*>
NanobenchmarkRegistry<Value_, Argument_>::NanobenchmarksMatching(
    std::regex const& filter) {
  std::vector<Nanobenchmark<Value_, Argument_> const*> matching;
  for (auto const& [name, nanobenchmark] :
       singleton().nanobenchmarks_by_name_) {
    if (std::regex_match(name, filter)) {
      matching.push_back(nanobenchmark);
    }
  }
  return matching;
}

template<typename Value_, typename Argument_>
Nanobenchmark<Value_, Argument_> const*
NanobenchmarkRegistry<Value_, Argument_>::NanobenchmarkFor(
    BenchmarkedFunction const function) {
  if (auto const it = singleton().nanobenchmarks_by_function_.find(function);
      it == singleton().nanobenchmarks_by_function_.end()) {
    return nullptr;
  } else {
    return it->second;
  }
}

template<typename Value_, typename Argument_>
NanobenchmarkRegistry<Value_, Argument_>&
NanobenchmarkRegistry<Value_, Argument_>::singleton() {
  static NanobenchmarkRegistry* instance = new NanobenchmarkRegistry;
  return *instance;
}

}  // namespace internal
}  // namespace _nanobenchmark
}  // namespace nanobenchmarks
}  // namespace principia
