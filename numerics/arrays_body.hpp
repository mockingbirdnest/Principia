#pragma once

#include "numerics/arrays.hpp"

namespace principia {
namespace numerics {

std::array<double, 2> const x{{1,2}};

template<typename Scalar, int size>
FixedVector<Scalar, size>::FixedVector() {
  Scalar zero{};  // Zero initialization.
  data_.fill(zero);
}

template<typename Scalar, int size>
FixedVector<Scalar, size>::FixedVector(std::array<Scalar, size> data)
    : data_(data) {}

template<typename Scalar, int size>
bool FixedVector<Scalar, size>::operator==(FixedVector const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows, int columns>
FixedMatrix<Scalar, rows, columns>::FixedMatrix(
    std::array<Scalar, rows * columns> data)
    : data_(data) {}

template<typename Scalar, int rows, int columns>
bool FixedMatrix<Scalar, rows, columns>::operator==(
    FixedMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows, int columns>
FixedVector<Scalar, rows> operator*(
    FixedMatrix<Scalar, rows, columns> const& left,
    FixedVector<Scalar, columns> const& right) {
  FixedVector<Scalar, rows> result;
  auto const* row_start = left.data_.data();
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      result.data_[i] += row_start[j] * right.data_[j];
    }
    row_start += columns;
  }
  return result;
}

}  // namespace numerics
}  // namespace principia
