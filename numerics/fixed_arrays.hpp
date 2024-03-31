#pragma once

#include <array>
#include <utility>
#include <vector>

#include "base/tags.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _fixed_arrays {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

template<typename Scalar_, int rows_, int columns_>
class FixedMatrix;
template<typename Scalar_, int rows_>
class FixedUpperTriangularMatrix;

template<typename Scalar_, int size_>
class FixedVector final {
 public:
  using Scalar = Scalar_;
  static constexpr int dimension = size_;

  constexpr FixedVector();
  explicit FixedVector(uninitialized_t);

  constexpr FixedVector(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedVector(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  // Convertible to an array.
  explicit constexpr operator std::array<Scalar, size_>() const;

  friend bool operator==(FixedVector const& left,
                         FixedVector const& right) = default;
  friend bool operator!=(FixedVector const& left,
                         FixedVector const& right) = default;

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  constexpr FixedVector& operator=(Scalar const (&right)[size_]);

  constexpr FixedVector& operator+=(FixedVector const& right);
  constexpr FixedVector& operator-=(FixedVector const& right);
  constexpr FixedVector& operator*=(double right);
  constexpr FixedVector& operator/=(double right);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  FixedVector<double, size_> Normalize() const;

  static constexpr int size() { return size_; }

  typename std::array<Scalar, size_>::const_iterator begin() const;
  typename std::array<Scalar, size_>::const_iterator end() const;

  template<typename H>
  friend H AbslHashValue(H h, FixedVector const& vector) {
    for (int index = 0; index < size_; ++index) {
      h = H::combine(std::move(h), vector.data_[index] / si::Unit<Scalar>);
    }
    return h;
  }

 private:
  std::array<Scalar, size_> data_;

  template<typename L, typename R, int s>
  friend constexpr Product<L, R> operator*(
      L* const left,
      FixedVector<R, s> const& right);
  template<typename L, typename R, int s>
  friend constexpr Product<L, R> operator*(
      TransposedView<FixedVector<L, s>> const& left,
      FixedVector<R, s> const& right);
  template<typename L, typename R, int r, int c>
  friend constexpr FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
};

template<typename Scalar_, int rows_, int columns_>
class FixedMatrix final {
  static constexpr int size_ = rows_ * columns_;

 public:
  using Scalar = Scalar_;
  static constexpr int rows() { return rows_; }
  static constexpr int columns() { return columns_; }

  constexpr FixedMatrix();
  explicit FixedMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedMatrix(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedMatrix(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  constexpr explicit FixedMatrix(
      TransposedView<FixedMatrix<Scalar, columns_, rows_>> const& view);

  friend bool operator==(FixedMatrix const& left,
                         FixedMatrix const& right) = default;
  friend bool operator!=(FixedMatrix const& left,
                         FixedMatrix const& right) = default;

  // For  0 < i < rows and 0 < j < columns, the entry a_ij is accessed as
  // |a(i, j)|.  if i and j do not satisfy these conditions, the expression
  // |a(i, j)| implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  constexpr FixedMatrix& operator=(Scalar const (&right)[size_]);

  constexpr FixedMatrix& operator+=(FixedMatrix const& right);
  constexpr FixedMatrix& operator-=(FixedMatrix const& right);
  constexpr FixedMatrix& operator*=(double right);
  constexpr FixedMatrix& operator/=(double right);

  constexpr FixedMatrix& operator*=(
      FixedMatrix<double, rows_, columns_> const& right)
    requires(rows_ == columns_);

  template<int r>
  Scalar const* row() const;

  Scalar FrobeniusNorm() const;

  // Applies the matrix as a bilinear form.  Present for compatibility with
  // |SymmetricBilinearForm|.  Prefer to use |TransposedView| and |operator*|.
  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>>
      operator()(FixedVector<LScalar, columns_> const& left,
                 FixedVector<RScalar, rows_> const& right) const;

  static FixedMatrix Identity();

 private:
  std::array<Scalar, size_> data_;

  template<typename L, typename R, int r, int c>
  friend constexpr FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
};

template<typename Scalar_, int rows_>
class FixedStrictlyLowerTriangularMatrix final {
  static constexpr int size_ = rows_ * (rows_ - 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr int rows() { return rows_; }
  static constexpr int columns() { return rows_; }

  constexpr FixedStrictlyLowerTriangularMatrix();
  explicit FixedStrictlyLowerTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedStrictlyLowerTriangularMatrix(
      std::array<Scalar, size_> const& data);

  friend bool operator==(FixedStrictlyLowerTriangularMatrix const& left,
                         FixedStrictlyLowerTriangularMatrix const& right) =
      default;
  friend bool operator!=(FixedStrictlyLowerTriangularMatrix const& left,
                         FixedStrictlyLowerTriangularMatrix const& right) =
      default;

  // For  0 ≤ j < i < rows, the entry a_ij is accessed as |a(i, j)|.
  // if i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  constexpr FixedStrictlyLowerTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

  template<int r>
  Scalar const* row() const;

 private:
  std::array<Scalar, size_> data_;
};

template<typename Scalar_, int rows_>
class FixedLowerTriangularMatrix final {
  static constexpr int size_ = rows_ * (rows_ + 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr int rows() { return rows_; }
  static constexpr int columns() { return rows_; }

  constexpr FixedLowerTriangularMatrix();
  explicit FixedLowerTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedLowerTriangularMatrix(
      std::array<Scalar, size_> const& data);

  explicit FixedLowerTriangularMatrix(
      TransposedView<FixedUpperTriangularMatrix<Scalar, rows_>> const& view);

  friend bool operator==(FixedLowerTriangularMatrix const& left,
                         FixedLowerTriangularMatrix const& right) = default;
  friend bool operator!=(FixedLowerTriangularMatrix const& left,
                         FixedLowerTriangularMatrix const& right) = default;

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as |a(i, j)|.
  // if i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  constexpr FixedLowerTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  std::array<Scalar, size_> data_;
};

template<typename Scalar_, int columns_>
class FixedUpperTriangularMatrix final {
  static constexpr int size_ = columns_ * (columns_ + 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr int rows() { return columns_; }
  static constexpr int columns() { return columns_; }

  constexpr FixedUpperTriangularMatrix();
  explicit FixedUpperTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedUpperTriangularMatrix(
      std::array<Scalar, size_> const& data);

  explicit FixedUpperTriangularMatrix(
      TransposedView<FixedLowerTriangularMatrix<Scalar, columns_>> const& view);

  friend bool operator==(FixedUpperTriangularMatrix const& left,
                         FixedUpperTriangularMatrix const& right) = default;
  friend bool operator!=(FixedUpperTriangularMatrix const& left,
                         FixedUpperTriangularMatrix const& right) = default;

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as |a(i, j)|.
  // if i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  constexpr FixedUpperTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This transposes a trapezoidal slice to make it column-major.
  static std::array<Scalar, size_> Transpose(
      std::array<Scalar, size_> const& data);

  std::array<Scalar, size_> data_;
};

// Prefer using the operator* that takes a TransposedView.
template<typename LScalar, typename RScalar, int size>
constexpr Product<LScalar, RScalar> InnerProduct(
    FixedVector<LScalar, size> const& left,
    FixedVector<RScalar, size> const& right);

template<typename Scalar, int size>
constexpr FixedVector<double, size> Normalize(
    FixedVector<Scalar, size> const& vector);

template<typename LScalar, typename RScalar, int size>
constexpr FixedMatrix<Product<LScalar, RScalar>, size, size> SymmetricProduct(
    FixedVector<LScalar, size> const& left,
    FixedVector<RScalar, size> const& right);

template<typename Scalar, int size>
constexpr FixedMatrix<Square<Scalar>, size, size> SymmetricSquare(
    FixedVector<Scalar, size> const& vector);

// Additive groups.

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size> operator+(
    FixedVector<Scalar, size> const& right);

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns> operator+(
    FixedMatrix<Scalar, rows, columns> const& right);

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size> operator-(
    FixedVector<Scalar, size> const& right);

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns> operator-(
    FixedMatrix<Scalar, rows, columns> const& right);

template<typename LScalar, typename RScalar, int size>
constexpr FixedVector<Sum<LScalar, RScalar>, size> operator+(
    FixedVector<LScalar, size> const& left,
    FixedVector<RScalar, size> const& right);

template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedMatrix<Sum<LScalar, RScalar>, rows, columns> operator+(
    FixedMatrix<LScalar, rows, columns> const& left,
    FixedMatrix<RScalar, rows, columns> const& right);

template<typename LScalar, typename RScalar, int size>
constexpr FixedVector<Difference<LScalar, RScalar>, size> operator-(
    FixedVector<LScalar, size> const& left,
    FixedVector<RScalar, size> const& right);

template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedMatrix<Difference<LScalar, RScalar>, rows, columns> operator-(
    FixedMatrix<LScalar, rows, columns> const& left,
    FixedMatrix<RScalar, rows, columns> const& right);

// Vector spaces.

template<typename LScalar, typename RScalar, int size>
constexpr FixedVector<Product<LScalar, RScalar>, size> operator*(
    LScalar const& left,
    FixedVector<RScalar, size> const& right);

template<typename LScalar, typename RScalar, int size>
constexpr FixedVector<Product<LScalar, RScalar>, size> operator*(
    FixedVector<LScalar, size> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns>
operator*(LScalar const& left,
          FixedMatrix<RScalar, rows, columns> const& right);

template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns>
operator*(FixedMatrix<LScalar, rows, columns> const& left,
          RScalar const& right);

template<typename LScalar, typename RScalar, int size>
constexpr FixedVector<Quotient<LScalar, RScalar>, size> operator/(
    FixedVector<LScalar, size> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedMatrix<Quotient<LScalar, RScalar>, rows, columns>
operator/(FixedMatrix<LScalar, rows, columns> const& left,
          RScalar const& right);

// Hilbert space and algebra.

// TODO(phl): We should have a RowView.
template<typename LScalar, typename RScalar, int size>
constexpr Product<LScalar, RScalar> operator*(
    LScalar* const left,
    FixedVector<RScalar, size> const& right);

template<typename LScalar, typename RScalar, int size>
constexpr Product<LScalar, RScalar> operator*(
    TransposedView<FixedVector<LScalar, size>> const& left,
    FixedVector<RScalar, size> const& right);

template<typename LScalar, typename RScalar, int lsize, int rsize>
constexpr FixedMatrix<Product<LScalar, RScalar>, lsize, rsize> operator*(
    FixedVector<LScalar, lsize> const& left,
    TransposedView<FixedVector<RScalar, rsize>> const& right);

template<typename LScalar, typename RScalar,
         int rows, int dimension, int columns>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns>
operator*(FixedMatrix<LScalar, rows, dimension> const& left,
          FixedMatrix<RScalar, dimension, columns> const& right);

template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedVector<Product<LScalar, RScalar>, rows> operator*(
    FixedMatrix<LScalar, rows, columns> const& left,
    FixedVector<RScalar, columns> const& right);

// Use this operator to multiply a row vector with a matrix.  We don't have an
// operator returning a TransposedView as that would cause dangling references.
template<typename LScalar, typename RScalar, int rows, int columns>
constexpr FixedVector<Product<LScalar, RScalar>, columns> operator*(
    TransposedView<FixedMatrix<LScalar, rows, columns>> const& left,
    FixedVector<RScalar, rows> const& right);

// Ouput.

template<typename Scalar, int size>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size> const& vector);

template<typename Scalar, int rows, int columns>
std::ostream& operator<<(std::ostream& out,
                         FixedMatrix<Scalar, rows, columns> const& matrix);

template<typename Scalar, int rows>
std::ostream& operator<<(
    std::ostream& out,
    FixedStrictlyLowerTriangularMatrix<Scalar, rows> const& matrix);

template<typename Scalar, int rows>
std::ostream& operator<<(
    std::ostream& out,
    FixedLowerTriangularMatrix<Scalar, rows> const& matrix);

template<typename Scalar, int columns>
std::ostream& operator<<(
    std::ostream& out,
    FixedUpperTriangularMatrix<Scalar, columns> const& matrix);

}  // namespace internal

using internal::FixedLowerTriangularMatrix;
using internal::FixedMatrix;
using internal::FixedStrictlyLowerTriangularMatrix;
using internal::FixedUpperTriangularMatrix;
using internal::FixedVector;

}  // namespace _fixed_arrays
}  // namespace numerics
}  // namespace principia

#include "numerics/fixed_arrays_body.hpp"
