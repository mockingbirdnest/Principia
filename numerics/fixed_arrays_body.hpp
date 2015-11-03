#pragma once

#include "numerics/fixed_arrays.hpp"

#include <algorithm>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace numerics {

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size>::FixedVector() {
  data_.fill({});  // Zero-initialized.
}

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size>::FixedVector(
    std::array<Scalar, size> const& data)
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

template<typename Scalar, int size>
Scalar& FixedVector<Scalar, size>::operator[](int const index) {
  return data_[index];
}

template<typename Scalar, int size>
constexpr Scalar const& FixedVector<Scalar, size>::operator[](
    int const index) const {
  return data_[index];
}

template<typename Scalar, int size>
FixedVector<Scalar, size>::operator std::vector<Scalar>() const {
  std::vector<Scalar> result(data_.size());
  std::copy(data_.begin(), data_.end(), result.begin());
  return result;
}

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>::FixedMatrix(
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

template<typename Scalar, int rows>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar, rows>::
    FixedStrictlyLowerTriangularMatrix(
        std::array<Scalar, kDimension> const& data)
    : data_(data) {}

template<typename Scalar, int rows>
FixedStrictlyLowerTriangularMatrix<Scalar, rows>::
    FixedStrictlyLowerTriangularMatrix(
        std::initializer_list<Scalar> const& data) {
  CHECK_EQ(kDimension, data.size());
  std::copy(data.begin(), data.end(), data_.begin());
}

template<typename Scalar, int rows>
bool FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator==(
    FixedStrictlyLowerTriangularMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows>
FixedStrictlyLowerTriangularMatrix<Scalar, rows>&
FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator=(
    std::initializer_list<Scalar> const& right) {
  CHECK_EQ(kDimension, right.size());
  std::copy(right.begin(), right.end(), data_.begin());
  return *this;
}

template<typename Scalar, int rows>
Scalar* FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator[](
    int const index) {
  return &data_[index * (index - 1) / 2];
}

template<typename Scalar, int rows>
constexpr Scalar const*
FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator[](
    int const index) const {
  return &data_[index * (index - 1) / 2];
}

template<typename Scalar, int rows>
int constexpr FixedStrictlyLowerTriangularMatrix<Scalar, rows>::kDimension;

}  // namespace numerics
}  // namespace principia
