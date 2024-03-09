#pragma once

#include "numerics/concepts.hpp"

namespace principia {
namespace numerics {
namespace _transposed_view {
namespace internal {

using namespace principia::numerics::_concepts;

template<typename T>
struct TransposedView {
  T const& transpose;

  constexpr int rows() const requires two_dimensional<T>;
  constexpr int columns() const requires two_dimensional<T>;
  constexpr int size() const requires one_dimensional<T>;

  constexpr typename T::Scalar& operator[](int index)
    requires one_dimensional<T>;
  constexpr typename T::Scalar const& operator[](int index) const
    requires one_dimensional<T>;

  constexpr typename T::Scalar& operator()(int row, int column)
    requires two_dimensional<T>;
  constexpr typename T::Scalar const& operator()(int row, int column) const
    requires two_dimensional<T>;
};

template<class T>
TransposedView(T) -> TransposedView<T>;

}  // namespace internal

using internal::TransposedView;

}  // namespace _transposed_view
}  // namespace numerics
}  // namespace principia

#include "numerics/transposed_view_body.hpp"
