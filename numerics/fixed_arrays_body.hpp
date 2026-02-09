#pragma once

#include "numerics/fixed_arrays.hpp"

#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

#include "glog/logging.h"
#include "numerics/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _fixed_arrays {
namespace internal {

using namespace principia::numerics::_elementary_functions;

// A helper class to compute the dot product of two arrays.  `LScalar` and
// `RScalar` are the types of the elements of the arrays.  `Left` and `Right`
// are the (deduced) types of the arrays.  They must both have an operator[].
// The third argument must be `std::make_index_sequence<size>`.
template<typename LScalar, typename RScalar, typename>
struct DotProduct;

template<typename LScalar, typename RScalar, std::int64_t... i>
struct DotProduct<LScalar, RScalar, std::index_sequence<i...>> {
  template<typename Left, typename Right>
  static Product<LScalar, RScalar> Compute(Left const& left,
                                           Right const& right);
};


template<typename LScalar, typename RScalar, std::int64_t... i>
template<typename Left, typename Right>
Product<LScalar, RScalar>
DotProduct<LScalar, RScalar, std::index_sequence<i...>>::Compute(
    Left const& left,
    Right const& right) {
  return ((left[i] * right[i]) + ...);
}

// The `data_` member is aggregate-initialized with an empty list initializer,
// which performs value initialization on the components.  For quantities this
// calls the default constructor, for non-class types this does
// zero-initialization.
template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>::FixedVector()
    : data_{} {
  if constexpr (use_heap) {
    data_ = std::make_unique<std::array<Scalar, size_>>({});
  }
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
FixedVector<Scalar_, size_, use_heap>::FixedVector(uninitialized_t) {
  if constexpr (use_heap) {
    data_ = std::make_unique<std::array<Scalar, size_>>();
  }
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>::FixedVector(
    std::array<Scalar, size_> const& data)
    : data_(MakeData(data)) {}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>::FixedVector(
    std::array<Scalar, size_>&& data)
    : data_(MakeData(std::move(data))) {}

template<typename Scalar_, std::int64_t size_, bool use_heap>
template<typename T>
  requires std::same_as<typename T::Scalar, Scalar_>
constexpr FixedVector<Scalar_, size_, use_heap>::FixedVector(
    ColumnView<T> const& view) : FixedVector(uninitialized) {
  CONSTEXPR_DCHECK(view.size() == size_);
  for (std::int64_t i = 0; i < size_; ++i) {
    (*this)[i] = view[i];
  }
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>::
operator std::array<Scalar_, size_>() const {
  return data();
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr Scalar_& FixedVector<Scalar_, size_, use_heap>::operator[](
    std::int64_t const index) {
  CONSTEXPR_DCHECK(0 <= index);
  CONSTEXPR_DCHECK(index < size());
  return data()[index];
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr Scalar_ const& FixedVector<Scalar_, size_, use_heap>::operator[](
    std::int64_t const index) const {
  CONSTEXPR_DCHECK(0 <= index);
  CONSTEXPR_DCHECK(index < size());
  return data()[index];
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>&
FixedVector<Scalar_, size_, use_heap>::operator=(Scalar const (&right)[size_]) {
  std::copy(right, right + size_, data().data());
  return *this;
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>&
FixedVector<Scalar_, size_, use_heap>::operator+=(FixedVector const& right) {
  for (std::int64_t i = 0; i < size(); ++i) {
    data()[i] += right.data()[i];
  }
  return *this;
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>&
FixedVector<Scalar_, size_, use_heap>::operator-=(FixedVector const& right) {
  for (std::int64_t i = 0; i < size(); ++i) {
    data()[i] -= right.data()[i];
  }
  return *this;
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>&
FixedVector<Scalar_, size_, use_heap>::operator*=(double const right) {
  for (auto& d : data()) {
    d *= right;
  }
  return *this;
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr FixedVector<Scalar_, size_, use_heap>&
FixedVector<Scalar_, size_, use_heap>::operator/=(double const right) {
  for (auto& d : data()) {
    d /= right;
  }
  return *this;
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
Scalar_ FixedVector<Scalar_, size_, use_heap>::Norm() const {
  return Sqrt(Norm²());
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
Square<Scalar_> FixedVector<Scalar_, size_, use_heap>::Norm²() const {
  return DotProduct<Scalar, Scalar, std::make_index_sequence<size_>>::Compute(
      data(), data());
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
template<bool uh>
FixedVector<double, size_, uh> FixedVector<Scalar_, size_, use_heap>::
Normalize() const {
  return *this / Norm();
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
typename std::array<Scalar_, size_>::const_iterator
FixedVector<Scalar_, size_, use_heap>::begin() const {
  return data().cbegin();
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
typename std::array<Scalar_, size_>::const_iterator
FixedVector<Scalar_, size_, use_heap>::end() const {
  return data().cend();
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
std::array<Scalar_, size_> const&
FixedVector<Scalar_, size_, use_heap>::data() const {
  if constexpr (use_heap) {
    return *data_;
  } else {
    return data_;
  }
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
std::array<Scalar_, size_>&
FixedVector<Scalar_, size_, use_heap>::data() {
  if constexpr (use_heap) {
    return *data_;
  } else {
    return data_;
  }
}

template<typename Scalar_, std::int64_t size_, bool use_heap>
constexpr typename FixedVector<Scalar_, size_, use_heap>::Data
FixedVector<Scalar_, size_, use_heap>::MakeData(auto&&... args) {
  if constexpr (use_heap) {
    return std::make_unique<std::array<Scalar_, size_>>(
        std::forward<decltype(args)>(args)...);
  } else {
    return std::array<Scalar, size_>{std::forward<decltype(args)>(args)...};
  }
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>::FixedMatrix()
    : data_{} {}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
FixedMatrix<Scalar_, rows_, columns_, use_heap>::FixedMatrix(uninitialized_t) {}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>::FixedMatrix(
    std::array<Scalar, size_> const& data)
    : data_(data) {}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>::FixedMatrix(
    std::array<Scalar, size_>&& data)
    : data_(std::move(data)) {}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
template<bool uh>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>::FixedMatrix(
    TransposedView<FixedMatrix<Scalar, columns_, rows_, uh>> const& view)
    : FixedMatrix(uninitialized) {
  for (std::int64_t i = 0; i < rows_; ++i) {
    for (std::int64_t j = 0; j < columns_; ++j) {
      (*this)(i, j) = view(i, j);
    }
  }
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr Scalar_& FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator()(
    std::int64_t const row, std::int64_t const column) {
  CONSTEXPR_DCHECK(0 <= row);
  CONSTEXPR_DCHECK(row < rows());
  CONSTEXPR_DCHECK(0 <= column);
  CONSTEXPR_DCHECK(column < columns());
  return data_[row * columns() + column];
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr Scalar_ const& FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator()(
    std::int64_t const row, std::int64_t const column) const {
  CONSTEXPR_DCHECK(0 <= row);
  CONSTEXPR_DCHECK(row < rows());
  CONSTEXPR_DCHECK(0 <= column);
  CONSTEXPR_DCHECK(column < columns());
  return data_[row * columns() + column];
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>&
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator=(
    Scalar const (&right)[size_]) {
  std::copy(right, right + size_, data_.data());
  return *this;
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_, bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>&
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator+=(FixedMatrix const& right) {
  for (std::int64_t i = 0; i < size_; ++i) {
    data_[i] += right.data_[i];
  }
  return *this;
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>&
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator-=(FixedMatrix const& right) {
  for (std::int64_t i = 0; i < size_; ++i) {
    data_[i] -= right.data_[i];
  }
  return *this;
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>&
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator*=(double const right) {
  for (auto& d : data_) {
    d *= right;
  }
  return *this;
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>&
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator/=(double const right) {
  for (auto& d : data_) {
    d /= right;
  }
  return *this;
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
constexpr FixedMatrix<Scalar_, rows_, columns_, use_heap>&
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator*=(
    FixedMatrix<double, rows_, columns_, use_heap> const& right)
  requires(rows_ == columns_) {
  // TODO(egg): We don't need to copy the whole matrix.
  return *this = *this * right;
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
template<std::int64_t r>
Scalar_ const* FixedMatrix<Scalar_, rows_, columns_, use_heap>::row() const {
  static_assert(r < rows_);
  return &data_[r * columns()];
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
Scalar_ FixedMatrix<Scalar_, rows_, columns_, use_heap>::FrobeniusNorm() const {
  Square<Scalar> Σᵢⱼaᵢⱼ²{};
  for (std::int64_t i = 0; i < rows(); ++i) {
    for (std::int64_t j = 0; j < columns(); ++j) {
      Σᵢⱼaᵢⱼ² += Pow<2>((*this)(i, j));
    }
  }
  return Sqrt(Σᵢⱼaᵢⱼ²);
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
template<typename LScalar, typename RScalar, bool luh, bool ruh>
Product<Scalar_, Product<LScalar, RScalar>>
FixedMatrix<Scalar_, rows_, columns_, use_heap>::operator()(
    FixedVector<LScalar, columns_, luh> const& left,
    FixedVector<RScalar, rows_, ruh> const& right) const {
  return TransposedView{left} * (*this * right);  // NOLINT
}

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_,
         bool use_heap>
FixedMatrix<Scalar_, rows_, columns_, use_heap>
FixedMatrix<Scalar_, rows_, columns_, use_heap>::Identity()
  requires(std::is_arithmetic_v<Scalar_> && rows_ == columns_) {
  FixedMatrix<Scalar, rows(), columns(), use_heap> m;
  for (std::int64_t i = 0; i < rows(); ++i) {
    m(i, i) = 1;
  }
  return m;
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::
    FixedStrictlyLowerTriangularMatrix()
    : data_{} {}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::
    FixedStrictlyLowerTriangularMatrix(uninitialized_t) {}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::
FixedStrictlyLowerTriangularMatrix(std::array<Scalar, size_> const& data)
    : data_(data) {}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
template<bool uh>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::
operator FixedMatrix<Scalar_, rows_, rows_, uh>() const {
  FixedMatrix<Scalar, rows_, rows_, uh> result;  // Initialized.
  for (std::int64_t i = 0; i < rows_; ++i) {
    for (std::int64_t j = 0; j < i; ++j) {
      result(i, j) = (*this)(i, j);
    }
  }
  return result;
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr Scalar_&
FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) {
  CONSTEXPR_DCHECK(0 <= column);
  CONSTEXPR_DCHECK(column < row);
  CONSTEXPR_DCHECK(row < rows());
  return data_[row * (row - 1) / 2 + column];
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr Scalar_ const&
FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) const {
  CONSTEXPR_DCHECK(0 <= column);
  CONSTEXPR_DCHECK(column < row);
  CONSTEXPR_DCHECK(row < rows());
  return data_[row * (row - 1) / 2 + column];
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>&
FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::operator=(
    Scalar const (&right)[size_]) {
  std::copy(right, right + size_, data_.data());
  return *this;
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
template<std::int64_t r>
Scalar_ const*
FixedStrictlyLowerTriangularMatrix<Scalar_, rows_, use_heap>::row() const {
  static_assert(r < rows_);
  return &data_[r * (r - 1) / 2];
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
FixedLowerTriangularMatrix()
    : data_{} {}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
FixedLowerTriangularMatrix(uninitialized_t) {}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
FixedLowerTriangularMatrix(std::array<Scalar, size_> const& data)
    : data_(data) {}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
template<bool uh>
FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
FixedLowerTriangularMatrix(
    TransposedView<
        FixedUpperTriangularMatrix<Scalar, rows_, uh>> const& view)
    : FixedLowerTriangularMatrix(uninitialized) {
  for (std::int64_t i = 0; i < rows(); ++i) {
    for (std::int64_t j = 0; j <= i; ++j) {
      (*this)(i, j) = view(i, j);
    }
  }
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
template<bool uh>
constexpr FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
operator FixedMatrix<Scalar_, rows_, rows_, uh>() const {
  FixedMatrix<Scalar, rows_, rows_, uh> result;  // Initialized.
  for (std::int64_t i = 0; i < rows_; ++i) {
    for (std::int64_t j = 0; j <= i; ++j) {
      result(i, j) = (*this)(i, j);
    }
  }
  return result;
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr Scalar_& FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) {
  CONSTEXPR_DCHECK(0 <= column);
  CONSTEXPR_DCHECK(column <= row);
  CONSTEXPR_DCHECK(row < rows());
  return data_[row * (row + 1) / 2 + column];
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr Scalar_ const& FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) const {
  CONSTEXPR_DCHECK(0 <= column);
  CONSTEXPR_DCHECK(column <= row);
  CONSTEXPR_DCHECK(row < rows());
  return data_[row * (row + 1) / 2 + column];
}

template<typename Scalar_, std::int64_t rows_, bool use_heap>
constexpr FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>&
FixedLowerTriangularMatrix<Scalar_, rows_, use_heap>::operator=(
    Scalar const (&right)[size_]) {
  std::copy(right, right + size_, data_.data());
  return *this;
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedStrictlyUpperTriangularMatrix()
    : data_{} {}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedStrictlyUpperTriangularMatrix(uninitialized_t) {}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedStrictlyUpperTriangularMatrix(std::array<Scalar, size_> const& data)
    : data_(Transpose(data)) {}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
template<bool uh>
FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedStrictlyUpperTriangularMatrix(
    TransposedView<
        FixedStrictlyLowerTriangularMatrix<Scalar, columns_, uh>> const& view)
    : FixedStrictlyUpperTriangularMatrix(uninitialized) {
  for (std::int64_t i = 0; i < rows(); ++i) {
    for (std::int64_t j = i; j < columns(); ++j) {
      (*this)(i, j) = view(i, j);
    }
  }
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
template<bool uh>
constexpr FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
operator FixedMatrix<Scalar_, columns_, columns_, uh>() const {
  FixedMatrix<Scalar, columns_, columns_, uh> result;  // Initialized.
  for (std::int64_t j = 0; j < columns_; ++j) {
    for (std::int64_t i = 0; i < j; ++i) {
      result(i, j) = (*this)(i, j);
    }
  }
  return result;
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr Scalar_&
FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) {
  CONSTEXPR_DCHECK(0 <= row);
  CONSTEXPR_DCHECK(row < column);
  CONSTEXPR_DCHECK(column < columns());
  return data_[column * (column - 1) / 2 + row];
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr Scalar_ const&
FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) const {
  CONSTEXPR_DCHECK(0 <= row);
  CONSTEXPR_DCHECK(row < column);
  CONSTEXPR_DCHECK(column < columns());
  return data_[column * (column - 1) / 2 + row];
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>&
FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::operator=(
    Scalar const (&right)[size_]) {
  std::copy(right, right + size_, data_.data());
  data_ = Transpose(data_);
  return *this;
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
auto FixedStrictlyUpperTriangularMatrix<Scalar_, columns_, use_heap>::Transpose(
    std::array<Scalar, size_> const& data)
    -> std::array<Scalar, size_> {
  std::array<Scalar, rows() * columns()> full;
  std::int64_t index = 0;
  for (std::int64_t row = 0; row < rows(); ++row) {
    for (std::int64_t column = row + 1; column < columns(); ++column) {
      full[row * columns() + column] = data[index];
      ++index;
    }
  }

  std::array<Scalar, size_> result;
  index = 0;
  for (std::int64_t column = 0; column < columns(); ++column) {
    for (std::int64_t row = 0; row < column; ++row) {
      result[index] = full[row * columns() + column];
      ++index;
    }
  }
  return result;
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedUpperTriangularMatrix()
    : data_{} {}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedUpperTriangularMatrix(uninitialized_t) {}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedUpperTriangularMatrix(std::array<Scalar, size_> const& data)
    : data_(Transpose(data)) {}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
template<bool uh>
FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
FixedUpperTriangularMatrix(
    TransposedView<
        FixedLowerTriangularMatrix<Scalar, columns_, uh>> const& view)
    : FixedUpperTriangularMatrix(uninitialized) {
  for (std::int64_t i = 0; i < rows(); ++i) {
    for (std::int64_t j = i; j < columns(); ++j) {
      (*this)(i, j) = view(i, j);
    }
  }
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
template<bool uh>
constexpr FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
operator FixedMatrix<Scalar_, columns_, columns_, uh>() const {
  FixedMatrix<Scalar, columns_, columns_, uh> result;  // Initialized.
  for (std::int64_t j = 0; j < columns_; ++j) {
    for (std::int64_t i = 0; i <= j; ++i) {
      result(i, j) = (*this)(i, j);
    }
  }
  return result;
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr Scalar_& FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) {
  CONSTEXPR_DCHECK(0 <= row);
  CONSTEXPR_DCHECK(row <= column);
  CONSTEXPR_DCHECK(column < columns());
  return data_[column * (column + 1) / 2 + row];
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr Scalar_ const&
FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::
operator()(std::int64_t const row, std::int64_t const column) const {
  CONSTEXPR_DCHECK(0 <= row);
  CONSTEXPR_DCHECK(row <= column);
  CONSTEXPR_DCHECK(column < columns());
  return data_[column * (column + 1) / 2 + row];
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
constexpr FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>&
FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::operator=(
    Scalar const (&right)[size_]) {
  std::copy(right, right + size_, data_.data());
  data_ = Transpose(data_);
  return *this;
}

template<typename Scalar_, std::int64_t columns_, bool use_heap>
auto FixedUpperTriangularMatrix<Scalar_, columns_, use_heap>::Transpose(
    std::array<Scalar, size_> const& data)
    -> std::array<Scalar, size_> {
  std::array<Scalar, rows() * columns()> full;
  std::int64_t index = 0;
  for (std::int64_t row = 0; row < rows(); ++row) {
    for (std::int64_t column = row; column < columns(); ++column) {
      full[row * columns() + column] = data[index];
      ++index;
    }
  }

  std::array<Scalar, size_> result;
  index = 0;
  for (std::int64_t column = 0; column < columns(); ++column) {
    for (std::int64_t row = 0; row <= column; ++row) {
      result[index] = full[row * columns() + column];
      ++index;
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool luh, bool ruh>
constexpr Product<LScalar, RScalar> InnerProduct(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right) {
  return DotProduct<LScalar, RScalar, std::make_index_sequence<size>>::Compute(
      left, right);
}

template<typename Scalar, std::int64_t size, bool uh, bool vuh>
constexpr FixedVector<double, size, uh> Normalize(
    FixedVector<Scalar, size, vuh> const& vector) {
  return vector / vector.Norm();
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool luh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, size, size, uh>
SymmetricProduct(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right) {
  FixedMatrix<Product<LScalar, RScalar>, size, size, uh> result(uninitialized);
  for (std::int64_t i = 0; i < size; ++i) {
    for (std::int64_t j = 0; j < i; ++j) {
      auto const r = 0.5 * (left[i] * right[j] + left[j] * right[i]);
      result(i, j) = r;
      result(j, i) = r;
    }
    result(i, i) = left[i] * right[i];
  }
  return result;
}

template<typename Scalar, std::int64_t size,
         bool uh, bool vuh>
constexpr FixedMatrix<Square<Scalar>, size, size, uh> SymmetricSquare(
    FixedVector<Scalar, size, vuh> const& vector) {
  FixedMatrix<Square<Scalar>, size, size, uh> result(uninitialized);
  for (std::int64_t i = 0; i < size; ++i) {
    for (std::int64_t j = 0; j < i; ++j) {
      auto const r =  vector[i] * vector[j];
      result(i, j) = r;
      result(j, i) = r;
    }
    result(i, i) = Pow<2>(vector[i]);
  }
  return result;
}

template<typename Scalar, std::int64_t size,
         bool uh, bool ruh>
constexpr FixedVector<Scalar, size, uh> operator+(
    FixedVector<Scalar, size, ruh> const& right) {
  if constexpr (uh == ruh) {
    return right;
  } else {
    std::array<Scalar, size> result_array;
    for (std::int64_t i = 0; i < size; ++i) {
      result_array[i] = right[i];
    }
    return FixedVector<Scalar, size, uh>(std::move(result_array));
  }
}

template<typename Scalar, std::int64_t rows, std::int64_t columns,
         bool uh, bool ruh>
constexpr FixedMatrix<Scalar, rows, columns, uh> operator+(
    FixedMatrix<Scalar, rows, columns, ruh> const& right) {
  if constexpr (uh == ruh) {
    return right;
  } else {
    FixedMatrix<Scalar, rows, columns, uh> result(uninitialized);
    for (std::int64_t i = 0; i < rows; ++i) {
      for (std::int64_t j = 0; j < columns; ++j) {
        result(i, j) = right(i, j);
      }
    }
    return result;
  }
}

template<typename Scalar, std::int64_t size,
         bool uh, bool ruh>
constexpr FixedVector<Scalar, size, uh> operator-(
    FixedVector<Scalar, size, ruh> const& right) {
  std::array<Scalar, size> result_array;
  for (std::int64_t i = 0; i < size; ++i) {
    result_array[i] = -right[i];
  }
  return FixedVector<Scalar, size, uh>(std::move(result_array));
}

template<typename Scalar, std::int64_t rows, std::int64_t columns,
         bool uh, bool ruh>
constexpr FixedMatrix<Scalar, rows, columns, uh> operator-(
    FixedMatrix<Scalar, rows, columns, ruh> const& right) {
  FixedMatrix<Scalar, rows, columns, uh> result(uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = -right(i, j);
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool luh, bool ruh>
constexpr FixedVector<Sum<LScalar, RScalar>, size, uh> operator+(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right) {
  std::array<Sum<LScalar, RScalar>, size> result;
  for (std::int64_t i = 0; i < size; ++i) {
    result[i] = left[i] + right[i];
  }
  return FixedVector<Sum<LScalar, RScalar>, size, uh>(std::move(result));
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool luh, bool ruh>
constexpr FixedMatrix<Sum<LScalar, RScalar>, rows, columns, uh> operator+(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    FixedMatrix<RScalar, rows, columns, ruh> const& right) {
  FixedMatrix<Sum<LScalar, RScalar>, rows, columns, uh> result(uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = left(i, j) + right(i, j);
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool luh, bool ruh>
constexpr FixedVector<Difference<LScalar, RScalar>, size, uh> operator-(
    FixedVector<LScalar, size, luh> const& left,
    FixedVector<RScalar, size, ruh> const& right) {
  std::array<Difference<LScalar, RScalar>, size> result;
  for (std::int64_t i = 0; i < size; ++i) {
    result[i] = left[i] - right[i];
  }
  return FixedVector<Difference<LScalar, RScalar>, size, uh>(std::move(result));
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool luh, bool ruh>
constexpr FixedMatrix<Difference<LScalar, RScalar>, rows, columns, uh> operator-(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    FixedMatrix<RScalar, rows, columns, ruh> const& right) {
  FixedMatrix<Difference<LScalar, RScalar>, rows, columns, uh> result(
      uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = left(i, j) - right(i, j);
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool ruh>
constexpr FixedVector<Product<LScalar, RScalar>, size, uh> operator*(
    LScalar const& left,
    FixedVector<RScalar, size, ruh> const& right) {
  std::array<Product<LScalar, RScalar>, size> result;
  for (std::int64_t i = 0; i < size; ++i) {
    result[i] = left * right[i];
  }
  return FixedVector<Product<LScalar, RScalar>, size, uh>(std::move(result));
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool luh>
constexpr FixedVector<Product<LScalar, RScalar>, size, uh> operator*(
    FixedVector<LScalar, size, luh> const& left,
    RScalar const& right) {
  std::array<Product<LScalar, RScalar>, size> result;
  for (std::int64_t i = 0; i < size; ++i) {
    result[i] = left[i] * right;
  }
  return FixedVector<Product<LScalar, RScalar>, size, uh>(std::move(result));
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh> operator*(
    LScalar const& left,
    FixedMatrix<RScalar, rows, columns, ruh> const& right) {
  FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh> result(
      uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = left * right(i, j);
    }
  }
  return result;
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool luh>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh> operator*(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    RScalar const& right) {
  FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh> result(
      uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = left(i, j) * right;
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool uh, bool luh>
constexpr FixedVector<Quotient<LScalar, RScalar>, size, uh> operator/(
    FixedVector<LScalar, size, luh> const& left,
    RScalar const& right) {
  FixedVector<Quotient<LScalar, RScalar>, size, uh> result(uninitialized);
  for (std::int64_t i = 0; i < left.size(); ++i) {
    result[i] = left[i] / right;
  }
  return result;
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool luh>
constexpr FixedMatrix<Quotient<LScalar, RScalar>, rows, columns, uh> operator/(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    RScalar const& right) {
  FixedMatrix<Quotient<LScalar, RScalar>, rows, columns, uh> result(
      uninitialized);
  for (std::int64_t i = 0; i < rows; ++i) {
    for (std::int64_t j = 0; j < columns; ++j) {
      result(i, j) = left(i, j) / right;
    }
  }
  return result;
}

template<typename LScalar, typename RScalar, std::int64_t size, bool ruh>
constexpr Product<LScalar, RScalar> operator*(
    LScalar* const left,
    FixedVector<RScalar, size, ruh> const& right) {
  return DotProduct<LScalar, RScalar, std::make_index_sequence<size>>::Compute(
      left, right.data_);
}

template<typename LScalar, typename RScalar, std::int64_t size,
         bool luh, bool ruh>
constexpr Product<LScalar, RScalar> operator*(
    TransposedView<FixedVector<LScalar, size, luh>> const& left,
    FixedVector<RScalar, size, ruh> const& right) {
  return DotProduct<LScalar, RScalar, std::make_index_sequence<size>>::Compute(
      left.transpose.data_, right.data_);
}

template<typename LScalar, typename RScalar,
         std::int64_t lsize, std::int64_t rsize,
         bool uh, bool luh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, lsize, rsize, uh> operator*(
    FixedVector<LScalar, lsize, luh> const& left,
    TransposedView<FixedVector<RScalar, rsize, ruh>> const& right) {
  FixedMatrix<Product<LScalar, RScalar>, lsize, rsize, uh> result(uninitialized);
  for (std::int64_t i = 0; i < lsize; ++i) {
    for (std::int64_t j = 0; j < rsize; ++j) {
      result(i, j) = left[i] * right[j];
    }
  }
  return result;
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t dimension, std::int64_t columns,
         bool uh, bool luh, bool ruh>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh>
operator*(FixedMatrix<LScalar, rows, dimension, luh> const& left,
          FixedMatrix<RScalar, dimension, columns, ruh> const& right) {
  FixedMatrix<Product<LScalar, RScalar>, rows, columns, uh> result{};
  for (std::int64_t i = 0; i < left.rows(); ++i) {
    for (std::int64_t j = 0; j < right.columns(); ++j) {
      for (std::int64_t k = 0; k < left.columns(); ++k) {
        result(i, j) += left(i, k) * right(k, j);
      }
    }
  }
  return result;
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool luh, bool ruh>
constexpr FixedVector<Product<LScalar, RScalar>, rows, uh> operator*(
    FixedMatrix<LScalar, rows, columns, luh> const& left,
    FixedVector<RScalar, columns, ruh> const& right) {
  std::array<Product<LScalar, RScalar>, rows> result;
  auto const* row = left.data_.data();
  for (std::int64_t i = 0; i < rows; ++i) {
    result[i] =
        DotProduct<LScalar, RScalar, std::make_index_sequence<columns>>::Compute(
            row, right.data_);
    row += columns;
  }
  return FixedVector<Product<LScalar, RScalar>, rows, uh>(std::move(result));
}

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns,
         bool uh, bool luh, bool ruh>
constexpr FixedVector<Product<LScalar, RScalar>, columns, uh> operator*(
    TransposedView<FixedMatrix<LScalar, rows, columns, luh>> const& left,
    FixedVector<RScalar, rows, ruh> const& right) {
  std::array<Product<LScalar, RScalar>, columns> result{};
  for (std::int64_t i = 0; i < columns; ++i) {
    auto& result_i = result[i];
    for (std::int64_t j = 0; j < rows; ++j) {
      result_i += left(i, j) * right[j];
    }
  }
  return FixedVector<Product<LScalar, RScalar>, columns, uh>(std::move(result));
}

template<typename Scalar, std::int64_t size, bool uh>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size, uh> const& vector) {
  std::stringstream s;
  for (std::int64_t i = 0; i < size; ++i) {
    s << (i == 0 ? "{" : "") << vector[i]
      << (i == size - 1 ? "}" : ", ");
  }
  out << s.str();
  return out;
}

template<typename Scalar, std::int64_t rows, std::int64_t columns, bool uh>
std::ostream& operator<<(std::ostream& out,
                         FixedMatrix<Scalar, rows, columns, uh> const& matrix) {
  out << "rows: " << rows << " columns: " << columns << "\n";
  for (std::int64_t i = 0; i < rows; ++i) {
    out << "{";
    for (std::int64_t j = 0; j < columns; ++j) {
      out << matrix(i, j);
      if (j < columns - 1) {
        out << ", ";
      }
    }
    out << "}\n";
  }
  return out;
}

template<typename Scalar, std::int64_t rows, bool uh>
std::ostream& operator<<(
    std::ostream& out,
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, uh> const& matrix) {
  out << "rows: " << matrix.rows() << "\n";
  for (std::int64_t i = 0; i < matrix.rows(); ++i) {
    out << "{";
    for (std::int64_t j = 0; j < i; ++j) {
      out << matrix(i, j);
      if (j < i - 1) {
        out << ", ";
      }
    }
    out << "}\n";
  }
  return out;
}

template<typename Scalar, std::int64_t rows, bool uh>
std::ostream& operator<<(
    std::ostream& out,
    FixedLowerTriangularMatrix<Scalar, rows, uh> const& matrix) {
  out << "rows: " << matrix.rows() << "\n";
  for (std::int64_t i = 0; i < matrix.rows(); ++i) {
    out << "{";
    for (std::int64_t j = 0; j <= i; ++j) {
      out << matrix(i, j);
      if (j < i) {
        out << ", ";
      }
    }
    out << "}\n";
  }
  return out;
}

template<typename Scalar, std::int64_t columns, bool uh>
std::ostream& operator<<(
    std::ostream& out,
    FixedUpperTriangularMatrix<Scalar, columns, uh> const& matrix) {
  out << "columns: " << matrix.columns() << "\n";
  for (std::int64_t i = 0; i < matrix.columns(); ++i) {
    out << "{";
    for (std::int64_t j = i; j < matrix.columns(); ++j) {
      if (j > i) {
        out << ", ";
      }
      out << matrix(i, j);
    }
    out << "}\n";
  }
  return out;
}

}  // namespace internal
}  // namespace _fixed_arrays
}  // namespace numerics
}  // namespace principia
