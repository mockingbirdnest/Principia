#pragma once

#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/flags/declare.h"
#include "mathematica/logger.hpp"
#include "nanobenchmarks/latency_distribution_table.hpp"

ABSL_DECLARE_FLAG(std::size_t, loop_iterations);
ABSL_DECLARE_FLAG(std::size_t, samples);
ABSL_DECLARE_FLAG(std::string, benchmark_filter);
ABSL_DECLARE_FLAG(std::string, log_to_mathematica);
ABSL_DECLARE_FLAG(double, input);

namespace principia {
namespace nanobenchmarks {
namespace _nanobenchmark {
namespace internal {

using namespace principia::mathematica::_logger;
using namespace principia::nanobenchmarks::_latency_distribution_table;

template<typename Value_ = double, typename Argument_ = double>
class Nanobenchmark {
 public:
  using Value = Value_;
  using Argument = Argument_;
  using BenchmarkedFunction = Value(*)(Argument);

  Nanobenchmark() = default;
  virtual ~Nanobenchmark() = default;

  // We disable inlining on this function so that the overhead is independent of
  // the callsite, and so that we actually call the benchmarked function via a
  // function pointer, instead of inlining it.
  __declspec(noinline) LatencyDistributionTable Run(Logger* logger) const;

  BenchmarkedFunction function() const;
  std::string const& name() const;

 protected:
  virtual Argument PrepareArgument(double x) const = 0;
  virtual double ConsumeValue(Value value) const = 0;
  virtual Value NanobenchmarkCase(Argument argument) const = 0;

  void SetFunction(BenchmarkedFunction function);
  void SetName(std::string_view name);

 private:
  BenchmarkedFunction function_ = nullptr;
  std::string name_;
};


template<typename Nanobenchmark_>
class NanobenchmarkRegistry {
 public:
  using Nanobenchmark = Nanobenchmark_;
  using BenchmarkedFunction = typename Nanobenchmark::BenchmarkedFunction;

  static Nanobenchmark const*
  Register(Nanobenchmark const* nanobenchmark);

  static std::vector<Nanobenchmark const*>
  NanobenchmarksMatching(std::regex const& filter);

  static Nanobenchmark const*
  NanobenchmarkFor(BenchmarkedFunction function);

 private:
  NanobenchmarkRegistry() = default;
  static NanobenchmarkRegistry& singleton();
  absl::flat_hash_map<BenchmarkedFunction, Nanobenchmark const*>
      nanobenchmarks_by_function_;
  absl::btree_map<std::string, Nanobenchmark const*> nanobenchmarks_by_name_;
};

#define NANOBENCHMARK_REGISTERED_NAME(line, n) registered_##line##_##n

#define NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method) \
  BaseClass##_##Method##_Nanobenchmark

#define NANOBENCHMARK_CONCAT_NAME(Function) \
  FunctionFixture##_##Function##_Nanobenchmark

#define NANOBENCHMARK_DECLARE_FIXTURE(BaseClass, Method)                      \
  class BaseClass##_##Method##_Nanobenchmark : public BaseClass {             \
   public:                                                                    \
    BaseClass##_##Method##_Nanobenchmark() : BaseClass() {                    \
      SetName(#BaseClass "/" #Method);                                        \
    }                                                                         \
                                                                              \
   protected:                                                                 \
    BaseClass::Value NanobenchmarkCase(BaseClass::Argument x) const override; \
  };

#define NANOBENCHMARK_DECLARE(Function)                      \
  class FunctionFixture##_##Function##_Nanobenchmark         \
      : public Nanobenchmark<> {                             \
   public:                                                   \
    FunctionFixture##_##Function##_Nanobenchmark() {         \
      SetName(#Function);                                    \
    }                                                        \
                                                             \
   protected:                                                \
    double PrepareArgument(double const x) const override {  \
      return x;                                              \
    }                                                        \
    double ConsumeValue(double const value) const override { \
      return value;                                          \
    }                                                        \
    double NanobenchmarkCase(double x) const override;       \
  };

#define NANOBENCHMARK_DECLARE_FUNCTION2(line, Function)                     \
  class FunctionFixture##_##line##_Nanobenchmark : public Nanobenchmark<> { \
   public:                                                                  \
    FunctionFixture##_##line##_Nanobenchmark() {                            \
      SetFunction(&Function);                                               \
      SetName(#Function);                                                   \
    }                                                                       \
                                                                            \
   protected:                                                               \
    double PrepareArgument(double const x) const override {                 \
      return x;                                                             \
    }                                                                       \
    double ConsumeValue(double const value) const override {                \
      return value;                                                         \
    }                                                                       \
    double NanobenchmarkCase(double const x) const override {               \
      return Function(x);                                                   \
    }                                                                       \
  };
#define NANOBENCHMARK_DECLARE_FUNCTION(line, Function) \
  NANOBENCHMARK_DECLARE_FUNCTION2(line, Function)

#define NANOBENCHMARK_DECLARE_REGISTERED(line, TestName, BaseClass)       \
  static BaseClass const* NANOBENCHMARK_REGISTERED_NAME(line, TestName) = \
      (NanobenchmarkRegistry<BaseClass>::Register(new TestName))

#define NANOBENCHMARK_REGISTER_FIXTURE(line, BaseClass, Method) \
  NANOBENCHMARK_DECLARE_REGISTERED(                             \
      line, NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method), BaseClass)

#define NANOBENCHMARK_REGISTER(line, Function) \
  NANOBENCHMARK_DECLARE_REGISTERED(            \
      line, NANOBENCHMARK_CONCAT_NAME(Function), Nanobenchmark<>)

#define NANOBENCHMARK_REGISTER_FUNCTION(line) \
  NANOBENCHMARK_DECLARE_REGISTERED(           \
      line, NANOBENCHMARK_CONCAT_NAME(line), Nanobenchmark<>)

#define NANOBENCHMARK_FIXTURE(BaseClass, Method, ...)          \
  NANOBENCHMARK_DECLARE_FIXTURE(BaseClass, Method)             \
  NANOBENCHMARK_REGISTER_FIXTURE(__LINE__, BaseClass, Method); \
  BaseClass::Value NANOBENCHMARK_CONCAT_NAME_FIXTURE(          \
      BaseClass,                                               \
      Method)::NanobenchmarkCase(BaseClass::Argument const argument) const

#define NANOBENCHMARK(Function)                                  \
  NANOBENCHMARK_DECLARE(Function)                                \
  NANOBENCHMARK_REGISTER(__LINE__, Function);                    \
  double NANOBENCHMARK_CONCAT_NAME(Function)::NanobenchmarkCase( \
      double const x) const

#define NANOBENCHMARK_FUNCTION(Function)             \
  NANOBENCHMARK_DECLARE_FUNCTION(__LINE__, Function) \
  NANOBENCHMARK_REGISTER_FUNCTION(__LINE__);

#define NANOBENCHMARK_EXTERN_C_FUNCTION(f) \
  extern "C" double f(double);             \
  NANOBENCHMARK_FUNCTION(f)

}  // namespace internal

using internal::Nanobenchmark;
using internal::NanobenchmarkRegistry;

}  // namespace _nanobenchmark
}  // namespace nanobenchmarks
}  // namespace principia
