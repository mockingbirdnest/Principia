#pragma once

#include <functional>
#include <map>
#include <regex>
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

 private:
  std::string name_;
};

class Fixture : public Nanobenchmark {
 public:
  Fixture() = default;

  virtual LatencyDistributionTable Run() override;

 protected:
  virtual double NanobenchmarkCase(double x) = 0;
};

class NanobenchmarkRegistry {
 public:
  static Nanobenchmark* Register(Nanobenchmark* nanobenchmark);

  static std::vector<Nanobenchmark*> NanobenchmarksMatching(
      std::regex const& filter);

 private:
  NanobenchmarkRegistry() = default;
  static NanobenchmarkRegistry& singleton();
  std::map<std::string, Nanobenchmark*> nanobenchmarks_by_name_;
};

#define NANOBENCHMARK_REGISTERED_NAME(line, n) registered_##line##_##n

#define NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method) \
  BaseClass##_##Method##_Nanobenchmark

#define NANOBENCHMARK_CONCAT_NAME(Function) \
  FunctionFixture##_##Function##_Nanobenchmark

#define NANOBENCHMARK_DECLARE_FIXTURE(BaseClass, Method)          \
  class BaseClass##_##Method##_Nanobenchmark : public BaseClass { \
   public:                                                        \
    BaseClass##_##Method##_Nanobenchmark() : BaseClass() {        \
      SetName(#BaseClass "/" #Method);                            \
    }                                                             \
                                                                  \
   protected:                                                     \
    double NanobenchmarkCase(double x) override;                  \
  };

#define NANOBENCHMARK_DECLARE(Function)                                 \
  class FunctionFixture##_##Function##_Nanobenchmark : public Fixture { \
   public:                                                              \
    FunctionFixture##_##Function##_Nanobenchmark() : Fixture() {        \
      SetName(#Function);                                               \
    }                                                                   \
                                                                        \
   protected:                                                           \
    double NanobenchmarkCase(double x) override;                        \
  };

#define NANOBENCHMARK_DECLARE_FUNCTION2(line, Function)             \
  class FunctionFixture##_##line##_Nanobenchmark : public Fixture { \
   public:                                                          \
    FunctionFixture##_##line##_Nanobenchmark() : Fixture() {        \
      SetName(#Function);                                           \
    }                                                               \
                                                                    \
   protected:                                                       \
    double NanobenchmarkCase(double x) override {                   \
      return Function(x);                                           \
    }                                                               \
  };
#define NANOBENCHMARK_DECLARE_FUNCTION(line, Function) \
  NANOBENCHMARK_DECLARE_FUNCTION2(line, Function)

#define NANOBENCHMARK_DECLARE_REGISTERED(line, TestName)                \
  static Nanobenchmark* NANOBENCHMARK_REGISTERED_NAME(line, TestName) = \
      (NanobenchmarkRegistry::Register(new TestName))

#define NANOBENCHMARK_REGISTER_FIXTURE(line, BaseClass, Method) \
  NANOBENCHMARK_DECLARE_REGISTERED(                             \
      line, NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method))

#define NANOBENCHMARK_REGISTER(line, Function) \
  NANOBENCHMARK_DECLARE_REGISTERED(line, NANOBENCHMARK_CONCAT_NAME(Function))

#define NANOBENCHMARK_REGISTER_FUNCTION(line) \
  NANOBENCHMARK_DECLARE_REGISTERED(line, NANOBENCHMARK_CONCAT_NAME(line))

#define NANOBENCHMARK_FIXTURE(BaseClass, Method, ...)          \
  NANOBENCHMARK_DECLARE_FIXTURE(BaseClass, Method)             \
  NANOBENCHMARK_REGISTER_FIXTURE(__LINE__, BaseClass, Method); \
  double NANOBENCHMARK_CONCAT_NAME_FIXTURE(                    \
      BaseClass, Method)::NanobenchmarkCase(double const x)

#define NANOBENCHMARK(Function)               \
  NANOBENCHMARK_DECLARE(Function)             \
  NANOBENCHMARK_REGISTER(__LINE__, Function); \
  double NANOBENCHMARK_CONCAT_NAME(Function)::NanobenchmarkCase(double const x)

#define NANOBENCHMARK_FUNCTION(Function)             \
  NANOBENCHMARK_DECLARE_FUNCTION(__LINE__, Function) \
  NANOBENCHMARK_REGISTER_FUNCTION(__LINE__);

#define NANOBENCHMARK_EXTERN_C_FUNCTION(f) \
  extern "C" double f(double);             \
  NANOBENCHMARK_FUNCTION(f)

}  // namespace internal

using internal::BenchmarkedFunction;
using internal::Fixture;
using internal::Nanobenchmark;
using internal::NanobenchmarkRegistry;

}  // namespace _function_registry
}  // namespace nanobenchmarks
}  // namespace principia
