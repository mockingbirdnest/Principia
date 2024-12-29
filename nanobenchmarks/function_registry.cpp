#include "nanobenchmarks/function_registry.hpp"

#include "glog/logging.h"

namespace principia {
namespace nanobenchmarks {
namespace _function_registry {
namespace internal {

bool FunctionRegistry::Register(std::string_view name,
                                BenchmarkedFunction function) {
  CHECK(names_by_function_.emplace(function, name).second)
      << " Registering function " << function << " as " << name << ": "
      << "function already registered as " << names_by_function_[function];
  CHECK(functions_by_name_.emplace(name, function).second)
      << " Registering function " << function << " as " << name << ": "
      << " name already taken by " << functions_by_name_.find(name)->second;
  return true;
}

FunctionRegistry& FunctionRegistry::singleton = *new FunctionRegistry();

std::map<std::string, BenchmarkedFunction, std::less<>> const&
FunctionRegistry::functions_by_name() const {
  return functions_by_name_;
}

std::map<BenchmarkedFunction, std::string> const&
FunctionRegistry::names_by_function() const {
  return names_by_function_;
}

}  // namespace internal
}  // namespace _function_registry
}  // namespace nanobenchmarks
}  // namespace principia