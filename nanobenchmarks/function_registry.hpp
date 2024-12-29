#pragma once

#include <map>
#include <string>
#include <string_view>

namespace principia {
namespace nanobenchmarks {
namespace _function_registry {
namespace internal {

#define BENCHMARK_CALLING_CONVENTION
using BenchmarkedFunction = double(BENCHMARK_CALLING_CONVENTION*)(double);

class FunctionRegistry {
 public:
  bool Register(std::string_view name, BenchmarkedFunction function);
  static FunctionRegistry& singleton;
  std::map<std::string, BenchmarkedFunction, std::less<>> const&
  functions_by_name() const;
  std::map<BenchmarkedFunction, std::string> const& names_by_function() const;

 private:
  FunctionRegistry() = default;
  std::map<std::string, BenchmarkedFunction, std::less<>> functions_by_name_;
  std::map<BenchmarkedFunction, std::string> names_by_function_;
};

#define BENCHMARK_FUNCTION(f)                                             \
  static bool registered_##f =                                            \
      ::principia::nanobenchmarks::_function_registry::FunctionRegistry:: \
          singleton.Register(#f, &(f))

#define EXPAND(x) x


#define BENCHMARK_FUNCTION_WITH_NAME(name, ...) \
  BENCHMARK_FUNCTION_WITH_NAME_INTERNAL(__LINE__, name, __VA_ARGS__)
#define BENCHMARK_FUNCTION_WITH_NAME_INTERNAL(line, name, ...) \
  BENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, __VA_ARGS__)
#define BENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, ...)           \
  namespace {                                                             \
  static bool registered##line =                                          \
      ::principia::nanobenchmarks::_function_registry::FunctionRegistry:: \
          singleton.Register(name, &(__VA_ARGS__));                       \
  }

#define BENCHMARKED_FUNCTION(f)                     \
  double BENCHMARK_CALLING_CONVENTION f(double x); \
  BENCHMARK_FUNCTION(f);                           \
  double BENCHMARK_CALLING_CONVENTION f(double x)

#define BENCHMARK_EXTERN_C_FUNCTION(f)                      \
  extern "C" double BENCHMARK_CALLING_CONVENTION f(double); \
  BENCHMARK_FUNCTION(f)

}  // namespace internal

using internal::BenchmarkedFunction;
using internal::FunctionRegistry;

}  // namespace _function_registry
}  // namespace nanobenchmarks
}  // namespace principia