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

template<typename T>
concept unbounded = requires(T& t) {
  { t.size() } -> std::same_as<int>;
} || requires(T& t) {
  { t.rows() } -> std::same_as<int>;
  { t.columns() } -> std::same_as<int>;
};

}  // namespace internal

using internal::one_dimensional;
using internal::two_dimensional;
using internal::unbounded;

}  // namespace _concepts
}  // namespace numerics
}  // namespace principia
