#pragma once

#include <concepts>
#include <type_traits>

namespace principia {
namespace numerics {
namespace _concepts {
namespace internal {

template<typename T>
concept one_dimensional = requires(T& t, int index) {
  { t[index] } -> std::same_as<typename T::Scalar&>;
};

template<typename T>
concept two_dimensional = requires(T& t, int row, int column) {
  { t(row, column) } -> std::same_as<typename T::Scalar&>;
};

}  // namespace internal

using internal::one_dimensional;
using internal::two_dimensional;

}  // namespace _concepts
}  // namespace numerics
}  // namespace principia