#pragma once

#include "numerics/arrays.hpp"

#include <algorithm>

#include "glog/logging.h"

namespace principia {
namespace numerics {

template<typename Scalar, int size>
FixedVector<Scalar, size>::FixedVector() {
  Scalar zero{};
  data_.fill(zero);
}

template<typename Scalar, int size>
FixedVector<Scalar, size>::FixedVector(std::array<Scalar, size> const& data)
    : data_(data) {}

template<typename Scalar, int size>
FixedVector<Scalar, size>::FixedVector(
    std::initializer_list<Scalar> const& data) {
  CHECK_EQ(size, data.size());
  std::copy(data.begin(), data.end(), data_.begin());
}

template<typename Scalar, int size>
bool FixedVector<Scalar, size>::operator==(FixedVector const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int size>
FixedVector<Scalar, size>& FixedVector<Scalar, size>::operator=(
    std::initializer_list<Scalar> const& right) {
  CHECK_EQ(size, right.size());
  std::copy(right.begin(), right.end(), data_.begin());
  return *this;
}

template<typename Scalar, int rows, int columns>
FixedMatrix<Scalar, rows, columns>::FixedMatrix(
    std::array<Scalar, rows * columns> const& data)
    : data_(data) {}

template<typename Scalar, int rows, int columns>
FixedMatrix<Scalar, rows, columns>::FixedMatrix(
    std::initializer_list<Scalar> const& data) {
  CHECK_EQ(rows * columns, data.size());
  std::copy(data.begin(), data.end(), data_.begin());
}

template<typename Scalar, int rows, int columns>
bool FixedMatrix<Scalar, rows, columns>::operator==(
    FixedMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows, int columns>
FixedMatrix<Scalar, rows, columns>&
FixedMatrix<Scalar, rows, columns>::operator=(
    std::initializer_list<Scalar> const& right) {
  CHECK_EQ(rows * columns, right.size());
  std::copy(right.begin(), right.end(), data_.begin());
  return *this;
}

template<typename ScalarLeft, typename ScalarRight, int rows, int columns>
FixedVector<Product<ScalarLeft, ScalarRight>, rows> operator*(
    FixedMatrix<ScalarLeft, rows, columns> const& left,
    FixedVector<ScalarRight, columns> const& right) {
  FixedVector<Product<ScalarLeft, ScalarRight>, rows> result;
  auto const* row = left.data_.data();
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      result.data_[i] += row[j] * right.data_[j];
    }
    row += columns;
  }
  return result;
}

}  // namespace numerics
}  // namespace principia
