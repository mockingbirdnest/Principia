#pragma once

#include "numerics/transposed_view.hpp"

namespace principia {
namespace numerics {
namespace _transposed_view {
namespace internal {

template<typename T>
constexpr int TransposedView<T>::rows() const
  requires two_dimensional<T> {
  // Note the transposition.
  return transpose.columns();
}

template<typename T>
constexpr int TransposedView<T>::columns() const
  requires two_dimensional<T> {
  // Note the transposition.
  return transpose.rows();
}

template<typename T>
constexpr int TransposedView<T>::size() const
  requires one_dimensional<T> {
  return transpose.size();
}

template<typename T>
constexpr typename T::Scalar& TransposedView<T>::operator[](int const index)
  requires one_dimensional<T> {
  return transpose[index];
}

template<typename T>
constexpr typename T::Scalar const& TransposedView<T>::operator[](
    int const index) const
  requires one_dimensional<T> {
  return transpose[index];
}

template<typename T>
constexpr typename T::Scalar& TransposedView<T>::operator()(int const row,
                                                            int const column)
  requires two_dimensional<T> {
  // Note the transposition.
  return transpose(column, row);
}

template<typename T>
constexpr typename T::Scalar const& TransposedView<T>::operator()(
    int const row,
    int const column) const
  requires two_dimensional<T> {
  // Note the transposition.
  return transpose(column, row);
}

}  // namespace internal
}  // namespace _transposed_view
}  // namespace numerics
}  // namespace principia
