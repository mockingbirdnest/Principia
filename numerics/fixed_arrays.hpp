#pragma once

#include <array>
#include <utility>
#include <vector>

#include "base/algebra.hpp"
#include "base/tags.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _fixed_arrays {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::base::_tags;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_si;

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_, bool use_heap_ = false>
class FixedMatrix;
template<typename Scalar_, std::int64_t rows_, bool use_heap_ = false>
class FixedUpperTriangularMatrix;

template<typename Scalar_, std::int64_t size_, bool use_heap_ = false>
class FixedVector final {
 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t dimension = size_;

  constexpr FixedVector();
  explicit FixedVector(uninitialized_t);

  constexpr FixedVector(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedVector(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  // Constructs a fixed vector by copying data from the view.  Note that the
  // result is fixed even if the matrix being viewed is an UnboundedMatrix.
  template<typename T>
    requires std::same_as<typename T::Scalar, Scalar_>
  constexpr explicit FixedVector(ColumnView<T> const& view);

  // Convertible to an array.
  explicit constexpr operator std::array<Scalar, size_>() const;

  friend bool operator==(FixedVector const& left,
                         FixedVector const& right) = default;
  friend bool operator!=(FixedVector const& left,
                         FixedVector const& right) = default;

  constexpr Scalar& operator[](std::int64_t index);
  constexpr Scalar const& operator[](std::int64_t index) const;

  constexpr FixedVector& operator=(Scalar const (&right)[size_]);

  constexpr FixedVector& operator+=(FixedVector const& right);
  constexpr FixedVector& operator-=(FixedVector const& right);
  constexpr FixedVector& operator*=(double right);
  constexpr FixedVector& operator/=(double right);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  FixedVector<double, size_> Normalize() const;

  static constexpr std::int64_t size() { return size_; }

  typename std::array<Scalar, size_>::const_iterator begin() const;
  typename std::array<Scalar, size_>::const_iterator end() const;

  template<typename H>
  friend H AbslHashValue(H h, FixedVector const& vector) {
    for (std::int64_t index = 0; index < size_; ++index) {
      h = H::combine(std::move(h), vector.data_[index] / si::Unit<Scalar>);
    }
    return h;
  }

 private:
  std::array<Scalar, size_> data_;

  template<typename L, typename R, std::int64_t s, bool uh>
  friend constexpr Product<L, R> operator*(
      L* const left,
      FixedVector<R, s, uh> const& right);
  template<typename L, typename R, std::int64_t s, bool uh>
  friend constexpr Product<L, R> operator*(
      TransposedView<FixedVector<L, s, uh>> const& left,
      FixedVector<R, s, uh> const& right);
  template<typename L, typename R, std::int64_t r, std::int64_t c, bool uh>
  friend constexpr FixedVector<Product<L, R>, r, uh> operator*(
      FixedMatrix<L, r, c, uh> const& left,
      FixedVector<R, c, uh> const& right);
};

template<typename Scalar_, std::int64_t rows_, std::int64_t columns_, bool use_heap_ = false>
class FixedMatrix final {
  static constexpr std::int64_t size_ = rows_ * columns_;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return rows_; }
  static constexpr std::int64_t columns() { return columns_; }

  constexpr FixedMatrix();
  explicit FixedMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedMatrix(
      std::array<Scalar, size_> const& data);  // NOLINT(runtime/explicit)
  constexpr FixedMatrix(
      std::array<Scalar, size_>&& data);  // NOLINT(runtime/explicit)

  constexpr explicit FixedMatrix(
      TransposedView<FixedMatrix<Scalar, columns_, rows_, use_heap_>> const& view);

  friend bool operator==(FixedMatrix const& left,
                         FixedMatrix const& right) = default;
  friend bool operator!=(FixedMatrix const& left,
                         FixedMatrix const& right) = default;

  // For  0 < i < rows and 0 < j < columns, the entry a_ij is accessed as
  // `a(i, j)`.  if i and j do not satisfy these conditions, the expression
  // `a(i, j)` implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  constexpr FixedMatrix& operator=(Scalar const (&right)[size_]);

  constexpr FixedMatrix& operator+=(FixedMatrix const& right);
  constexpr FixedMatrix& operator-=(FixedMatrix const& right);
  constexpr FixedMatrix& operator*=(double right);
  constexpr FixedMatrix& operator/=(double right);

  constexpr FixedMatrix& operator*=(
      FixedMatrix<double, rows_, columns_, use_heap_> const& right)
    requires(rows_ == columns_);

  template<std::int64_t r>
  Scalar const* row() const;

  Scalar FrobeniusNorm() const;

  // Applies the matrix as a bilinear form.  Present for compatibility with
  // `SymmetricBilinearForm`.  Prefer to use `TransposedView` and `operator*`.
  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>>
      operator()(FixedVector<LScalar, columns_, use_heap_> const& left,
                 FixedVector<RScalar, rows_, use_heap_> const& right) const;

  static FixedMatrix Identity()
    requires(std::is_arithmetic_v<Scalar_> && rows_ == columns_);

 private:
  std::array<Scalar, size_> data_;

  template<typename L, typename R, std::int64_t r, std::int64_t c, bool uh>
  friend constexpr FixedVector<Product<L, R>, r, uh> operator*(
      FixedMatrix<L, r, c, uh> const& left,
      FixedVector<R, c, uh> const& right);
};

template<typename Scalar_, std::int64_t rows_, bool use_heap_ = false>
class FixedStrictlyLowerTriangularMatrix final {
  static constexpr std::int64_t size_ = rows_ * (rows_ - 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return rows_; }
  static constexpr std::int64_t columns() { return rows_; }

  constexpr FixedStrictlyLowerTriangularMatrix();
  explicit FixedStrictlyLowerTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedStrictlyLowerTriangularMatrix(
      std::array<Scalar, size_> const& data);

  explicit constexpr operator FixedMatrix<Scalar, rows_, rows_, use_heap_>() const;

  friend bool operator==(FixedStrictlyLowerTriangularMatrix const& left,
                         FixedStrictlyLowerTriangularMatrix const& right) =
      default;
  friend bool operator!=(FixedStrictlyLowerTriangularMatrix const& left,
                         FixedStrictlyLowerTriangularMatrix const& right) =
      default;

  // For  0 ≤ j < i < rows, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  constexpr FixedStrictlyLowerTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

  template<std::int64_t r>
  Scalar const* row() const;

 private:
  std::array<Scalar, size_> data_;
};

template<typename Scalar_, std::int64_t rows_, bool use_heap_ = false>
class FixedLowerTriangularMatrix final {
  static constexpr std::int64_t size_ = rows_ * (rows_ + 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return rows_; }
  static constexpr std::int64_t columns() { return rows_; }

  constexpr FixedLowerTriangularMatrix();
  explicit FixedLowerTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedLowerTriangularMatrix(
      std::array<Scalar, size_> const& data);

  explicit FixedLowerTriangularMatrix(
      TransposedView<FixedUpperTriangularMatrix<Scalar, rows_, use_heap_>> const& view);

  explicit constexpr operator FixedMatrix<Scalar, rows_, rows_, use_heap_>() const;

  friend bool operator==(FixedLowerTriangularMatrix const& left,
                         FixedLowerTriangularMatrix const& right) = default;
  friend bool operator!=(FixedLowerTriangularMatrix const& left,
                         FixedLowerTriangularMatrix const& right) = default;

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  constexpr FixedLowerTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  std::array<Scalar, size_> data_;
};

template<typename Scalar_, std::int64_t columns_, bool use_heap_ = false>
class FixedStrictlyUpperTriangularMatrix final {
  static constexpr std::int64_t size_ = columns_ * (columns_ - 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return columns_; }
  static constexpr std::int64_t columns() { return columns_; }

  constexpr FixedStrictlyUpperTriangularMatrix();
  explicit FixedStrictlyUpperTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedStrictlyUpperTriangularMatrix(
      std::array<Scalar, size_> const& data);

  explicit FixedStrictlyUpperTriangularMatrix(
      TransposedView<
          FixedStrictlyLowerTriangularMatrix<Scalar, columns_, use_heap_>> const& view);

  explicit constexpr operator FixedMatrix<Scalar, columns_, columns_, use_heap_>() const;

  friend bool operator==(
      FixedStrictlyUpperTriangularMatrix const& left,
      FixedStrictlyUpperTriangularMatrix const& right) = default;
  friend bool operator!=(
      FixedStrictlyUpperTriangularMatrix const& left,
      FixedStrictlyUpperTriangularMatrix const& right) = default;

  // For  0 ≤ i < j < columns, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

  constexpr FixedStrictlyUpperTriangularMatrix& operator=(
      Scalar const (&right)[size_]);

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This transposes a trapezoidal slice to make it column-major.
  static std::array<Scalar, size_> Transpose(
      std::array<Scalar, size_> const& data);

  std::array<Scalar, size_> data_;
};

template<typename Scalar_, std::int64_t columns_, bool use_heap_ = false>
class FixedUpperTriangularMatrix final {
  static constexpr std::int64_t size_ = columns_ * (columns_ + 1) / 2;

 public:
  using Scalar = Scalar_;
  static constexpr std::int64_t rows() { return columns_; }
  static constexpr std::int64_t columns() { return columns_; }

  constexpr FixedUpperTriangularMatrix();
  explicit FixedUpperTriangularMatrix(uninitialized_t);

  // The `data` must be in row-major format.
  constexpr FixedUpperTriangularMatrix(
      std::array<Scalar, size_> const& data);

  explicit FixedUpperTriangularMatrix(
      TransposedView<FixedLowerTriangularMatrix<Scalar, columns_, use_heap_>> const& view);

  explicit constexpr operator FixedMatrix<Scalar, columns_, columns_, use_heap_>() const;

  friend bool operator==(FixedUpperTriangularMatrix const& left,
                         FixedUpperTriangularMatrix const& right) = default;
  friend bool operator!=(FixedUpperTriangularMatrix const& left,
                         FixedUpperTriangularMatrix const& right) = default;

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as `a(i, j)`.
  // if i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  constexpr Scalar& operator()(std::int64_t row, std::int64_t column);
  constexpr Scalar const& operator()(std::int64_t row,
                                     std::int64_t column) const;

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
template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr Product<LScalar, RScalar> InnerProduct(
    FixedVector<LScalar, size, use_heap> const& left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename Scalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<double, size, use_heap> Normalize(
    FixedVector<Scalar, size, use_heap> const& vector);

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr FixedMatrix<Product<LScalar, RScalar>, size, size, use_heap> SymmetricProduct(
    FixedVector<LScalar, size, use_heap> const& left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename Scalar, std::int64_t size, bool use_heap = false>
constexpr FixedMatrix<Square<Scalar>, size, size, use_heap> SymmetricSquare(
    FixedVector<Scalar, size, use_heap> const& vector);

// Additive groups.

template<typename Scalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Scalar, size, use_heap> operator+(
    FixedVector<Scalar, size, use_heap> const& right);

template<typename Scalar, std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Scalar, rows, columns, use_heap> operator+(
    FixedMatrix<Scalar, rows, columns, use_heap> const& right);

template<typename Scalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Scalar, size, use_heap> operator-(
    FixedVector<Scalar, size, use_heap> const& right);

template<typename Scalar, std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Scalar, rows, columns, use_heap> operator-(
    FixedMatrix<Scalar, rows, columns, use_heap> const& right);

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Sum<LScalar, RScalar>, size, use_heap> operator+(
    FixedVector<LScalar, size, use_heap> const& left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Sum<LScalar, RScalar>, rows, columns, use_heap> operator+(
    FixedMatrix<LScalar, rows, columns, use_heap> const& left,
    FixedMatrix<RScalar, rows, columns, use_heap> const& right);

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Difference<LScalar, RScalar>, size, use_heap> operator-(
    FixedVector<LScalar, size, use_heap> const& left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Difference<LScalar, RScalar>, rows, columns, use_heap> operator-(
    FixedMatrix<LScalar, rows, columns, use_heap> const& left,
    FixedMatrix<RScalar, rows, columns, use_heap> const& right);

// Vector spaces.

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Product<LScalar, RScalar>, size, use_heap> operator*(
    LScalar const& left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Product<LScalar, RScalar>, size, use_heap> operator*(
    FixedVector<LScalar, size, use_heap> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, use_heap>
operator*(LScalar const& left,
          FixedMatrix<RScalar, rows, columns, use_heap> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, use_heap>
operator*(FixedMatrix<LScalar, rows, columns, use_heap> const& left,
          RScalar const& right);

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr FixedVector<Quotient<LScalar, RScalar>, size, use_heap> operator/(
    FixedVector<LScalar, size, use_heap> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Quotient<LScalar, RScalar>, rows, columns, use_heap>
operator/(FixedMatrix<LScalar, rows, columns, use_heap> const& left,
          RScalar const& right);

// Hilbert space and algebra.

// TODO(phl): We should have a RowView.
template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr Product<LScalar, RScalar> operator*(
    LScalar* const left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename LScalar, typename RScalar, std::int64_t size, bool use_heap = false>
constexpr Product<LScalar, RScalar> operator*(
    TransposedView<FixedVector<LScalar, size, use_heap>> const& left,
    FixedVector<RScalar, size, use_heap> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t lsize, std::int64_t rsize, bool use_heap = false>
constexpr FixedMatrix<Product<LScalar, RScalar>, lsize, rsize, use_heap> operator*(
    FixedVector<LScalar, lsize, use_heap> const& left,
    TransposedView<FixedVector<RScalar, rsize, use_heap>> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t dimension, std::int64_t columns, bool use_heap = false>
constexpr FixedMatrix<Product<LScalar, RScalar>, rows, columns, use_heap>
operator*(FixedMatrix<LScalar, rows, dimension, use_heap> const& left,
          FixedMatrix<RScalar, dimension, columns, use_heap> const& right);

template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedVector<Product<LScalar, RScalar>, rows, use_heap> operator*(
    FixedMatrix<LScalar, rows, columns, use_heap> const& left,
    FixedVector<RScalar, columns, use_heap> const& right);

// Use this operator to multiply a row vector with a matrix.  We don't have an
// operator returning a TransposedView as that would cause dangling references.
template<typename LScalar, typename RScalar,
         std::int64_t rows, std::int64_t columns, bool use_heap = false>
constexpr FixedVector<Product<LScalar, RScalar>, columns, use_heap> operator*(
    TransposedView<FixedMatrix<LScalar, rows, columns, use_heap>> const& left,
    FixedVector<RScalar, rows, use_heap> const& right);

// Ouput.

template<typename Scalar, std::int64_t size, bool use_heap = false>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size, use_heap> const& vector);

template<typename Scalar, std::int64_t rows, std::int64_t columns, bool use_heap = false>
std::ostream& operator<<(std::ostream& out,
                         FixedMatrix<Scalar, rows, columns, use_heap> const& matrix);

template<typename Scalar, std::int64_t rows, bool use_heap = false>
std::ostream& operator<<(
    std::ostream& out,
    FixedStrictlyLowerTriangularMatrix<Scalar, rows, use_heap> const& matrix);

template<typename Scalar, std::int64_t rows, bool use_heap = false>
std::ostream& operator<<(
    std::ostream& out,
    FixedLowerTriangularMatrix<Scalar, rows, use_heap> const& matrix);

template<typename Scalar, std::int64_t columns, bool use_heap = false>
std::ostream& operator<<(
    std::ostream& out,
    FixedUpperTriangularMatrix<Scalar, columns, use_heap> const& matrix);

}  // namespace internal

using internal::FixedLowerTriangularMatrix;
using internal::FixedMatrix;
using internal::FixedStrictlyLowerTriangularMatrix;
using internal::FixedStrictlyUpperTriangularMatrix;
using internal::FixedUpperTriangularMatrix;
using internal::FixedVector;

}  // namespace _fixed_arrays
}  // namespace numerics
}  // namespace principia

#include "numerics/fixed_arrays_body.hpp"
