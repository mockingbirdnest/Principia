#include "functions/multiprecision.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {
namespace internal {

cpp_bin_float_50 Sin(cpp_rational const& angle) {
  return static_cast<cpp_bin_float_50>(angle);
}

}  // namespace internal
}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
