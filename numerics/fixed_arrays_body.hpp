
#pragma once

#include "numerics/fixed_arrays.hpp"

#include <algorithm>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace numerics {
namespace internal_fixed_arrays {

// A helper class to compute the dot product of two arrays.  |ScalarLeft| and
// |ScalarRight| are the types of the elements of the arrays.  |Left| and
// |Right| are the (deduced) types of the arrays.  They must both have an
// operator[].  |size| is the size of the arrays.
template<typename ScalarLeft, typename ScalarRight, int size, int i = size - 1>
struct DotProduct {
  template<typename Left, typename Right>
  static Product<ScalarLeft, ScalarRight> Compute(Left const& left,
                                                  Right const& right);
};

template<typename ScalarLeft, typename ScalarRight, int size>
struct DotProduct<ScalarLeft, ScalarRight, size, 0> {
  template<typename Left, typename Right>
  static Product<ScalarLeft, ScalarRight> Compute(Left const& left,
                                                  Right const& right);
};

template<typename ScalarLeft, typename ScalarRight, int size, int i>
template<typename Left, typename Right>
Product<ScalarLeft, ScalarRight>
DotProduct<ScalarLeft, ScalarRight, size, i>::Compute(Left const& left,
                                                      Right const& right) {
  return left[i] * right[i] +
         DotProduct<ScalarLeft, ScalarRight, size, i - 1>::Compute(left, right);
}

template<typename ScalarLeft, typename ScalarRight, int size>
template<typename Left, typename Right>
Product<ScalarLeft, ScalarRight>
DotProduct<ScalarLeft, ScalarRight, size, 0>::Compute(Left const& left,
                                                      Right const& right) {
  return left[0] * right[0];
}

// The |data_| member is aggregate-initialized with an empty list initializer,
// which performs value initialization on the components.  For quantities this
// calls the default constructor, for non-class types this does
// zero-initialization.
template<typename Scalar, int size_>
constexpr FixedVector<Scalar, size_>::FixedVector() : data_{} {}

template<typename Scalar, int size_>
FixedVector<Scalar, size_>::FixedVector(uninitialized_t) {}

template<typename Scalar, int size_>
constexpr FixedVector<Scalar, size_>::FixedVector(
    std::array<Scalar, size_> const& data)
    : data_(data) {}

template<typename Scalar, int size_>
constexpr FixedVector<Scalar, size_>::FixedVector(
    std::array<Scalar, size_>&& data)
    : data_(std::move(data)) {}

template<typename Scalar, int size_>
bool FixedVector<Scalar, size_>::operator==(FixedVector const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int size_>
constexpr Scalar& FixedVector<Scalar, size_>::operator[](int const index) {
  return data_[index];
}

template<typename Scalar, int size_>
constexpr Scalar const& FixedVector<Scalar, size_>::operator[](
    int const index) const {
  return data_[index];
}

template<typename Scalar, int size_>
FixedVector<Scalar, size_>::operator std::vector<Scalar>() const {
  std::vector<Scalar> result(data_.size());
  std::copy(data_.begin(), data_.end(), result.begin());
  return result;
}

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>::FixedMatrix()
    : data_{} {}

template<typename Scalar, int rows, int columns>
FixedMatrix<Scalar, rows, columns>::FixedMatrix(uninitialized_t) {}

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>::FixedMatrix(
    std::array<Scalar, rows * columns> const& data)
    : data_(data) {}

template<typename Scalar, int rows, int columns>
bool FixedMatrix<Scalar, rows, columns>::operator==(
    FixedMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows, int columns>
Scalar* FixedMatrix<Scalar, rows, columns>::operator[](int const index) {
  return &data_[index * columns];
}

template<typename Scalar, int rows, int columns>
constexpr Scalar const* FixedMatrix<Scalar, rows, columns>::operator[](
    int const index) const {
  return &data_[index * columns];
}

template<typename Scalar, int rows, int columns>
template<int r>
FixedMatrix<Scalar, rows, columns>::Row<r>::Row(const FixedMatrix* const matrix)
    : matrix_(matrix) {}

template<typename Scalar, int rows, int columns>
template<int r>
constexpr Scalar const& FixedMatrix<Scalar, rows, columns>::Row<r>::operator[](
    int index) const {
  return (matrix_->data_)[r * columns + index];
}

template<typename Scalar, int rows, int columns>
template<int r>
template<typename S>
Product<Scalar, S>
FixedMatrix<Scalar, rows, columns>::Row<r>::operator*(
    FixedVector<S, columns> const& right) {
  return DotProduct<Scalar, S, columns>::Compute(*this, right);
}

template<typename Scalar, int rows, int columns>
template<int r>
typename FixedMatrix<Scalar, rows, columns>::template Row<r>
FixedMatrix<Scalar, rows, columns>::row() const {
  return Row<r>(this);
}

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr FixedVector<Difference<ScalarLeft, ScalarRight>, size> operator-(
    FixedVector<ScalarLeft, size> const& left,
    FixedVector<ScalarRight, size> const& right) {
  std::array<Difference<ScalarLeft, ScalarRight>, size> result{};
  for (int i = 0; i < size; ++i) {
    result[i] = left[i] - right[i];
  }
  return FixedVector<Difference<ScalarLeft, ScalarRight>, size>(
      std::move(result));
}

template<typename ScalarLeft, typename ScalarRight, int rows, int columns>
FixedVector<Product<ScalarLeft, ScalarRight>, rows> operator*(
    FixedMatrix<ScalarLeft, rows, columns> const& left,
    FixedVector<ScalarRight, columns> const& right) {
  std::array<Product<ScalarLeft, ScalarRight>, rows> result;
  auto const* row = left.data_.data();
  for (int i = 0; i < rows; ++i) {
    result[i] =
        DotProduct<ScalarLeft, ScalarRight, columns>::Compute(row, right.data_);
    row += columns;
  }
  return FixedVector<Product<ScalarLeft, ScalarRight>, rows>(std::move(result));
}

template<typename Scalar, int size>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size> const& vector) {
  std::stringstream s;
  for (int i = 0; i < size; ++i) {
    s << (i == 0 ? "{" : "") << vector.data_[i]
      << (i == size - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}

template<typename Scalar, int rows>
std::ostream& operator<<(
    std::ostream& out,
    FixedLowerTriangularMatrix<Scalar, rows> const& matrix) {
  out << "rows: " << matrix.rows << "\n";
  for (int i = 0; i < matrix.rows; ++i) {
    out << "{";
    for (int j = 0; j <= i; ++j) {
      out << matrix[i][j];
      if (j < i) {
        out << ", ";
      }
    }
    out << "}\n";
  }
  return out;
}

template<typename Scalar, int columns>
std::ostream& operator<<(
    std::ostream& out,
    FixedUpperTriangularMatrix<Scalar, columns> const& matrix) {
  out << "columns: " << matrix.columns() << "\n";
  for (int i = 0; i < matrix.columns(); ++i) {
    out << "{";
    for (int j = i; j < matrix.columns(); ++j) {
      if (j > i) {
        out << ", ";
      }
      out << matrix[i][j];
    }
    out << "}\n";
  }
  return out;
}

template<typename Scalar, int rows_>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::
    FixedStrictlyLowerTriangularMatrix()
    : data_{} {}

template<typename Scalar, int rows_>
FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::
    FixedStrictlyLowerTriangularMatrix(uninitialized_t) {}

template<typename Scalar, int rows_>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::
    FixedStrictlyLowerTriangularMatrix(
        std::array<Scalar, dimension> const& data)
    : data_(data) {}

template<typename Scalar, int rows_>
bool FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::operator==(
    FixedStrictlyLowerTriangularMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows_>
Scalar* FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::operator[](
    int const index) {
  return &data_[index * (index - 1) / 2];
}

template<typename Scalar, int rows_>
constexpr Scalar const*
FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::operator[](
    int const index) const {
  return &data_[index * (index - 1) / 2];
}

template<typename Scalar, int rows_>
constexpr int FixedStrictlyLowerTriangularMatrix<Scalar, rows_>::dimension;

template<typename Scalar, int rows_>
constexpr FixedLowerTriangularMatrix<Scalar, rows_>::
FixedLowerTriangularMatrix()
    : data_{} {}

template<typename Scalar, int rows_>
FixedLowerTriangularMatrix<Scalar, rows_>::FixedLowerTriangularMatrix(
    uninitialized_t) {}

template<typename Scalar, int rows_>
constexpr FixedLowerTriangularMatrix<Scalar, rows_>::
    FixedLowerTriangularMatrix(std::array<Scalar, dimension> const& data)
    : data_(data) {}

template<typename Scalar, int rows_>
bool FixedLowerTriangularMatrix<Scalar, rows_>::operator==(
    FixedLowerTriangularMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows_>
Scalar* FixedLowerTriangularMatrix<Scalar, rows_>::operator[](
    int const index) {
  return &data_[index * (index + 1) / 2];
}

template<typename Scalar, int rows_>
constexpr Scalar const*
FixedLowerTriangularMatrix<Scalar, rows_>::operator[](
    int const index) const {
  return &data_[index * (index + 1) / 2];
}

template<typename Scalar, int rows_>
constexpr int FixedLowerTriangularMatrix<Scalar, rows_>::dimension;

template<typename Scalar, int columns_>
constexpr FixedUpperTriangularMatrix<Scalar, columns_>::
FixedUpperTriangularMatrix()
    : data_{} {}

template<typename Scalar, int columns_>
FixedUpperTriangularMatrix<Scalar, columns_>::FixedUpperTriangularMatrix(
    uninitialized_t) {}

template<typename Scalar, int columns_>
constexpr FixedUpperTriangularMatrix<Scalar, columns_>::
    FixedUpperTriangularMatrix(std::array<Scalar, dimension> const& data)
    : data_(Transpose(data)) {}

template<typename Scalar, int columns_>
bool FixedUpperTriangularMatrix<Scalar, columns_>::operator==(
    FixedUpperTriangularMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int columns_>
template<typename Matrix>
Scalar& FixedUpperTriangularMatrix<Scalar, columns_>::Row<Matrix>::operator[](
    int const column) {
  DCHECK_LT(column, columns_);
  return matrix_.data_[column * (column + 1) / 2 + row_];
}

template<typename Scalar, int columns_>
template<typename Matrix>
Scalar const&
FixedUpperTriangularMatrix<Scalar, columns_>::Row<Matrix>::operator[](
    int const column) const {
  DCHECK_LT(column, columns_);
  return matrix_.data_[column * (column + 1) / 2 + row_];
}

template<typename Scalar, int columns_>
template<typename Matrix>
FixedUpperTriangularMatrix<Scalar, columns_>::Row<Matrix>::Row(Matrix& matrix,
                                                               int const row)
    : matrix_(const_cast<std::remove_const_t<Matrix>&>(matrix)),
      row_(row) {}

template<typename Scalar, int columns_>
auto FixedUpperTriangularMatrix<Scalar, columns_>::operator[](int const row)
    -> Row<FixedUpperTriangularMatrix<Scalar, columns_>> {
  return Row<FixedUpperTriangularMatrix<Scalar, columns_>>{*this, row};
}

template<typename Scalar, int columns_>
auto FixedUpperTriangularMatrix<Scalar, columns_>::operator[](int const row)
    const -> Row<FixedUpperTriangularMatrix<Scalar, columns_> const> {
  return Row<FixedUpperTriangularMatrix<Scalar, columns_> const>{*this, row};
}

template<typename Scalar, int columns_>
auto FixedUpperTriangularMatrix<Scalar, columns_>::Transpose(
    std::array<Scalar, dimension> const& data)
    -> std::array<Scalar, dimension> {
  std::array<Scalar, columns_ * columns_> full;
  int index = 0;
  for (int row = 0; row < columns_; ++row) {
    for (int column = row; column < columns_; ++column) {
      full[row * columns_ + column] = data[index];
      ++index;
    }
  }

  std::array<Scalar, dimension> result;
  index = 0;
  for (int column = 0; column < columns_; ++column) {
    for (int row = 0; row <= column; ++row) {
      result[index] = full[row * columns_ + column];
      ++index;
    }
  }
  return result;
}

template<typename Scalar, int columns_>
constexpr int FixedUpperTriangularMatrix<Scalar, columns_>::dimension;

}  // namespace internal_fixed_arrays
}  // namespace numerics
}  // namespace principia
