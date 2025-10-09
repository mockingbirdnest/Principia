#pragma once

#include <functional>
#include <map>
#include <string>
#include <string_view>

#include "absl/flags/declare.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_COMPILER_MSVC.
#include "nanobenchmarks/flag_parsing.hpp"
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

  std::string const& name() const;

 protected:
  void SetName(std::string_view name);

  std::string name_;
};

class Fixture : public Nanobenchmark {
 public:
  Fixture() = default;

  virtual LatencyDistributionTable Run() override;

 protected:
  virtual double NanobenchmarkCase(double x) = 0;
};

class FunctionNanobenchmark : public Fixture {
 public:
  FunctionNanobenchmark(std::string_view name, BenchmarkedFunction function);

  virtual LatencyDistributionTable Run() override;

 private:
  BenchmarkedFunction function_;
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

class NanobenchmarkRegistry {
 public:
  static Nanobenchmark* Register(Nanobenchmark* nanobenchmark);
  static std::map<std::string, Nanobenchmark*> const& nanobenchmarks_by_name();

 private:
  NanobenchmarkRegistry() = default;
  static NanobenchmarkRegistry& singleton();
  std::map<std::string, Nanobenchmark*> nanobenchmarks_by_name_;
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
      SetName(#BaseClass "/" #Method);                            \
    }                                                             \
                                                                  \
   protected:                                                     \
    virtual double NanobenchmarkCase(double x) override;          \
  };

#define NANOBENCHMARK_PRIVATE_REGISTER_F(line, TestName) \
  NANOBENCHMARK_PRIVATE_DECLARE(line, TestName) =        \
      (NanobenchmarkRegistry::Register(new TestName))

#define NANOBENCHMARK_REGISTER_F(line, BaseClass, Method) \
  NANOBENCHMARK_PRIVATE_REGISTER_F(                       \
      line, NANOBENCHMARK_PRIVATE_CONCAT_NAME(BaseClass, Method))

#define NANOBENCHMARK(n)                       \
  NANOBENCHMARK_PRIVATE_DECLARE(__LINE__, n) = \
      (NanobenchmarkRegistry::Register(new FunctionNanobenchmark(#n, n)))

#define NANOBENCHMARK_F(BaseClass, Method, ...)          \
  NANOBENCHMARK_PRIVATE_DECLARE_F(BaseClass, Method)     \
  NANOBENCHMARK_REGISTER_F(__LINE__, BaseClass, Method); \
  double NANOBENCHMARK_PRIVATE_CONCAT_NAME(              \
      BaseClass, Method)::NanobenchmarkCase(double const x)

////Obsolete

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
