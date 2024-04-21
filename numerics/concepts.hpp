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
} || requires(T const& t, int row, int column) {
  { t(row, column) } -> std::same_as<typename T::Scalar const&>;
};

template<typename T1, typename T2>
concept same_elements_as =
    std::same_as<typename T1::Scalar, typename T2::Scalar>;

}  // namespace internal

using internal::one_dimensional;
using internal::same_elements_as;
using internal::two_dimensional;

}  // namespace _concepts
}  // namespace numerics
}  // namespace principia
