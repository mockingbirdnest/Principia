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

// TODO(phl): This should support the same operations as unbounded_arrays.hpp.

template<typename Scalar, int rows, int columns>
class FixedMatrix;

template<typename Scalar, int size_>
class FixedVector final {
 public:
  static constexpr int dimension = size_;

  constexpr FixedVector();
  explicit FixedVector(uninitialized_t);

  // TODO(egg): Figure out why we have a move-conversion for |FixedVector| but
  // not the matrices.
  constexpr FixedVector(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedVector(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  static constexpr int size() { return size_; }

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  typename std::array<Scalar, size_>::const_iterator begin() const;
  typename std::array<Scalar, size_>::const_iterator end() const;

  bool operator==(FixedVector const& right) const;
  bool operator!=(FixedVector const& right) const;

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

template<typename Scalar, int rows_, int columns_>
class FixedMatrix final {
 public:
  static constexpr int rows() { return rows_; }
  static constexpr int columns() { return columns_; }
  static constexpr int size() { return rows_ * columns_; }

  constexpr FixedMatrix();
  explicit FixedMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedMatrix(
      std::array<Scalar, size()> const& data);  // NOLINT(runtime/explicit)

  // For  0 < i < rows and 0 < j < columns, the entry a_ij is accessed as
  // |a(i, j)|.  if i and j do not satisfy these conditions, the expression
  // |a(i, j)| implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  template<int r>
  Scalar const* row() const;

  FixedMatrix Transpose() const;

  Scalar FrobeniusNorm() const;

  bool operator==(FixedMatrix const& right) const;
  bool operator!=(FixedMatrix const& right) const;

  // Applies the matrix as a bilinear form.  Present for compatibility with
  // |SymmetricBilinearForm|.  Prefer to use |TransposedView| and |operator*|.
  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>>
      operator()(FixedVector<LScalar, columns_> const& left,
                 FixedVector<RScalar, rows_> const& right) const;

  static FixedMatrix Identity();

 private:
  std::array<Scalar, size()> data_;

  template<typename L, typename R, int r, int c>
  friend constexpr FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
};

template<typename Scalar, int rows_>
class FixedStrictlyLowerTriangularMatrix final {
 public:
  static constexpr int rows() { return rows_; }
  static constexpr int columns() { return rows_; }
  static constexpr int size() { return rows_ * (rows_ - 1) / 2; }

  constexpr FixedStrictlyLowerTriangularMatrix();
  explicit FixedStrictlyLowerTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedStrictlyLowerTriangularMatrix(
      std::array<Scalar, size()> const& data);

  // For  0 ≤ j < i < rows, the entry a_ij is accessed as |a(i, j)|.
  // if i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  template<int r>
  Scalar const* row() const;

  bool operator==(FixedStrictlyLowerTriangularMatrix const& right) const;
  bool operator!=(FixedStrictlyLowerTriangularMatrix const& right) const;

 private:
  std::array<Scalar, size()> data_;
};

template<typename Scalar, int rows_>
class FixedLowerTriangularMatrix final {
 public:
  static constexpr int rows() { return rows_; }
  static constexpr int columns() { return rows_; }
  static constexpr int size() { return rows_ * (rows_ + 1) / 2; }

  constexpr FixedLowerTriangularMatrix();
  explicit FixedLowerTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedLowerTriangularMatrix(
      std::array<Scalar, size()> const& data);

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as |a(i, j)|.
  // if i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  bool operator==(FixedLowerTriangularMatrix const& right) const;
  bool operator!=(FixedLowerTriangularMatrix const& right) const;

 private:
  std::array<Scalar, size()> data_;
};

template<typename Scalar, int columns_>
class FixedUpperTriangularMatrix final {
 public:
  static constexpr int rows() { return columns_; }
  static constexpr int columns() { return columns_; }
  static constexpr int size() { return columns_ * (columns_ + 1) / 2; }

  constexpr FixedUpperTriangularMatrix();
  explicit FixedUpperTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedUpperTriangularMatrix(
      std::array<Scalar, size()> const& data);

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as |a(i, j)|.
  // if i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  constexpr Scalar& operator()(int row, int column);
  constexpr Scalar const& operator()(int row, int column) const;

  bool operator==(FixedUpperTriangularMatrix const& right) const;
  bool operator!=(FixedUpperTriangularMatrix const& right) const;

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This transposes a trapezoidal slice to make it column-major.
  static std::array<Scalar, size()> Transpose(
      std::array<Scalar, size()> const& data);

  std::array<Scalar, size()> data_;
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

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size>& operator+=(
    FixedVector<Scalar, size>& left,
    FixedVector<Scalar, size> const& right);

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>& operator+=(
    FixedMatrix<Scalar, rows, columns>& left,
    FixedMatrix<Scalar, rows, columns> const& right);

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size>& operator-=(
    FixedVector<Scalar, size>& left,
    FixedVector<Scalar, size> const& right);

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>& operator-=(
    FixedMatrix<Scalar, rows, columns>& left,
    FixedMatrix<Scalar, rows, columns> const& right);

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

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size>& operator*=(
    FixedVector<Scalar, size>& left,
    double right);

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>& operator*=(
    FixedMatrix<Scalar, rows, columns>& left,
    double right);

template<typename Scalar, int size>
constexpr FixedVector<Scalar, size>& operator/=(
    FixedVector<Scalar, size>& left,
    double right);

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>& operator/=(
    FixedMatrix<Scalar, rows, columns>& left,
    double right);

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
