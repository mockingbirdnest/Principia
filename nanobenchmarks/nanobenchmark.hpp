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
  virtual Value NanobenchmarkCase(Argument argument) const = 0;

  void SetFunction(BenchmarkedFunction function);
  void SetName(std::string_view name);

 private:
  BenchmarkedFunction function_ = nullptr;
  std::string name_;
};


template<typename Value_, typename Argument_>
class NanobenchmarkRegistry {
 public:
  using BenchmarkedFunction =
      typename Nanobenchmark<Value_, Argument_>::BenchmarkedFunction;

  static Nanobenchmark<Value_, Argument_> const* Register(
      Nanobenchmark<Value_, Argument_> const* nanobenchmark);

  static std::vector<Nanobenchmark<Value_, Argument_> const*>
  NanobenchmarksMatching(std::regex const& filter);

  static Nanobenchmark<Value_, Argument_> const* NanobenchmarkFor(
      BenchmarkedFunction function);

 private:
  NanobenchmarkRegistry() = default;
  static NanobenchmarkRegistry& singleton();
  absl::flat_hash_map<BenchmarkedFunction,
                      Nanobenchmark<Value_, Argument_> const*>
      nanobenchmarks_by_function_;
  absl::btree_map<std::string, Nanobenchmark<Value_, Argument_> const*>
      nanobenchmarks_by_name_;
};

#define NANOBENCHMARK_REGISTERED_NAME(line, n) registered_##line##_##n

#define NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method) \
  BaseClass##_##Method##_Nanobenchmark

#define NANOBENCHMARK_CONCAT_NAME(Function) \
  FunctionFixture##_##Function##_Nanobenchmark

#define NANOBENCHMARK_DECLARE_FIXTURE(BaseClass, Method)                      \
  class NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method)                  \
      : public BaseClass {                                                    \
   public:                                                                    \
    NANOBENCHMARK_CONCAT_NAME_FIXTURE(BaseClass, Method)() : BaseClass() {    \
      SetName(#BaseClass "/" #Method);                                        \
    }                                                                         \
                                                                              \
   protected:                                                                 \
    BaseClass::Value NanobenchmarkCase(BaseClass::Argument x) const override; \
  };

#define NANOBENCHMARK_DECLARE(Function)                                \
  class NANOBENCHMARK_CONCAT_NAME(Function) : public Nanobenchmark<> { \
   public:                                                             \
    NANOBENCHMARK_CONCAT_NAME(Function)() {                            \
      SetName(#Function);                                              \
    }                                                                  \
                                                                       \
   protected:                                                          \
    Value NanobenchmarkCase(Argument x) const override;                \
  };

#define NANOBENCHMARK_DECLARE_FUNCTION2(line, Function)            \
  class NANOBENCHMARK_CONCAT_NAME(line) : public Nanobenchmark<> { \
   public:                                                         \
    NANOBENCHMARK_CONCAT_NAME(line)() {                            \
      SetFunction(&Function);                                      \
      SetName(#Function);                                          \
    }                                                              \
                                                                   \
   protected:                                                      \
    Value NanobenchmarkCase(Argument const x) const override {     \
      return Function(x);                                          \
    }                                                              \
  };
#define NANOBENCHMARK_DECLARE_FUNCTION(line, Function) \
  NANOBENCHMARK_DECLARE_FUNCTION2(line, Function)

#define NANOBENCHMARK_DECLARE_REGISTERED(line, FixtureClass)                \
  static auto const* NANOBENCHMARK_REGISTERED_NAME(line, FixtureClass) =    \
      (NanobenchmarkRegistry<FixtureClass::Value, FixtureClass::Argument>:: \
           Register(new FixtureClass))

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

#include "nanobenchmarks/nanobenchmark_body.hpp"
