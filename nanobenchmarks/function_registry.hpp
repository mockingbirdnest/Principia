#pragma once

#include <functional>
#include <map>
#include <string>
#include <string_view>

namespace principia {
namespace nanobenchmarks {
namespace _function_registry {
namespace internal {

using BenchmarkedFunction = double (*)(double);

class FunctionRegistry {
 public:
  static bool Register(std::string_view name, BenchmarkedFunction function);
  static std::map<std::string, BenchmarkedFunction, std::less<>> const&
  functions_by_name();
  static std::map<BenchmarkedFunction, std::string> const& names_by_function();

 private:
  FunctionRegistry() = default;
  static FunctionRegistry& singleton_;
  std::map<std::string, BenchmarkedFunction, std::less<>> functions_by_name_;
  std::map<BenchmarkedFunction, std::string> names_by_function_;
};

#define BENCHMARK_FUNCTION_WITH_NAME(name, ...) \
  BENCHMARK_FUNCTION_WITH_NAME_INTERNAL(__LINE__, name, __VA_ARGS__)
#define BENCHMARK_FUNCTION_WITH_NAME_INTERNAL(line, name, ...) \
  BENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, __VA_ARGS__)
#define BENCHMARK_FUNCTION_WITH_NAME_INTERNAL2(line, name, ...)             \
  namespace {                                                               \
  static bool registered##line = ::principia::nanobenchmarks::              \
      _function_registry::FunctionRegistry::Register(name, &(__VA_ARGS__)); \
  }


#define BENCHMARK_FUNCTION(...) \
  BENCHMARK_FUNCTION_WITH_NAME(#__VA_ARGS__, __VA_ARGS__)

#define BENCHMARKED_FUNCTION(f) \
  double f(double x);           \
  BENCHMARK_FUNCTION(f);        \
  double f(double x)

#define BENCHMARK_EXTERN_C_FUNCTION(f) \
  extern "C" double f(double);         \
  BENCHMARK_FUNCTION(f)

}  // namespace internal

using internal::BenchmarkedFunction;
using internal::FunctionRegistry;

}  // namespace _function_registry
}  // namespace nanobenchmarks
}  // namespace principia
