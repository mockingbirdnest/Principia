#pragma once

#include <functional>
#include <map>
#include <string>
#include <string_view>

#include "absl/flags/declare.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "nanobenchmarks/latency_distribution_table.hpp"

ABSL_DECLARE_FLAG(std::size_t, loop_iterations);
ABSL_DECLARE_FLAG(std::size_t, samples);
ABSL_DECLARE_FLAG(std::string, benchmark_filter);
ABSL_DECLARE_FLAG(std::string, log_to_mathematica);
ABSL_DECLARE_FLAG(double, input);
ABSL_DECLARE_FLAG(std::vector<double>, quantiles);

namespace principia {
namespace nanobenchmarks {
namespace _function_registry {
namespace internal {

using namespace principia::nanobenchmarks::_latency_distribution_table;

using BenchmarkedFunction = double (*)(double);

class Nanobenchmark {
 public:
  Nanobenchmark() = default;
  virtual ~Nanobenchmark() = default;
  virtual LatencyDistributionTable Run() = 0;

 protected:
  void SetName(std::string const& name) {
    name_ = name;
  }

  std::string name_;
};

class Fixture : public Nanobenchmark {
 public:
  Fixture() = default;
  virtual ~Fixture() = default;

  virtual LatencyDistributionTable Run() override {
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

 protected:
  virtual double NanobenchmarkCase(double x) = 0;
};

class FunctionRegistry {
 public:
  static bool Register(std::string_view name, BenchmarkedFunction function);
  static std::map<std::string, BenchmarkedFunction, std::less<>> const&
  functions_by_name();
  static std::map<BenchmarkedFunction, std::string> const& names_by_function();

 private:
  FunctionRegistry() = default;
  static FunctionRegistry& singleton();
  std::map<std::string, BenchmarkedFunction, std::less<>> functions_by_name_;
  std::map<BenchmarkedFunction, std::string> names_by_function_;
};

#define NANOBENCHMARK_FUNCTION_WITH_NAME(name, ...) \
  NANOBENCHMARK_FUNCTION_WITH_NAME_INTERNAL(__LINE__, name, __VA_ARGS__)
#define NANOBENCHMARK_FUNCTION_WITH_NAME_INTERNAL(line, name, ...) \
  NANOBENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, __VA_ARGS__)
#define NANOBENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, ...)         \
  namespace {                                                               \
  static bool registered##line = ::principia::nanobenchmarks::              \
      _function_registry::FunctionRegistry::Register(name, &(__VA_ARGS__)); \
  }

#define NANOBENCHMARK_PRIVATE_NAME(line, n) registered_##line##_##n

#define NANOBENCHMARK_PRIVATE_DECLARE(line, n) \
  static Nanobenchmark* NANOBENCHMARK_PRIVATE_NAME(line, n)

#define NANOBENCHMARK_PRIVATE_CONCAT_NAME(BaseClass, Method) \
  BaseClass##_##Method##_Nanobenchmark

#define NANOBENCHMARK_PRIVATE_DECLARE_F(BaseClass, Method)        \
  class BaseClass##_##Method##_Nanobenchmark : public BaseClass { \
   public:                                                        \
    BaseClass##_##Method##_Nanobenchmark() : BaseClass() {        \
      this->SetName(#BaseClass "/" #Method);                      \
    }                                                             \
                                                                  \
   protected:                                                     \
    virtual double NanobenchmarkCase(double x) override;          \
  };

#define NANOBENCHMARK_PRIVATE_REGISTER_F(line, TestName) \
  NANOBENCHMARK_PRIVATE_DECLARE(line, TestName) = nullptr;
//(RegisterNanobenchmarkInternal(new TestName))

#define NANOBENCHMARK_REGISTER_F(line, BaseClass, Method) \
  NANOBENCHMARK_PRIVATE_REGISTER_F(                       \
      line, NANOBENCHMARK_PRIVATE_CONCAT_NAME(BaseClass, Method))

#define NANOBENCHMARK_F(BaseClass, Method, ...)          \
  NANOBENCHMARK_PRIVATE_DECLARE_F(BaseClass, Method)     \
  NANOBENCHMARK_REGISTER_F(__LINE__, BaseClass, Method); \
  double NANOBENCHMARK_PRIVATE_CONCAT_NAME(              \
      BaseClass, Method)::NanobenchmarkCase(double const x)

#define NANOBENCHMARK_FUNCTION(...) \
  NANOBENCHMARK_FUNCTION_WITH_NAME(#__VA_ARGS__, __VA_ARGS__)

#define NANOBENCHMARKED_FUNCTION(f) \
  double f(double x);               \
  NANOBENCHMARK_FUNCTION(f);        \
  double f(double x)

#define NANOBENCHMARK_EXTERN_C_FUNCTION(f) \
  extern "C" double f(double);             \
  NANOBENCHMARK_FUNCTION(f)

}  // namespace internal

using internal::BenchmarkedFunction;
using internal::Fixture;
using internal::FunctionRegistry;
using internal::Nanobenchmark;

}  // namespace _function_registry
}  // namespace nanobenchmarks
}  // namespace principia
