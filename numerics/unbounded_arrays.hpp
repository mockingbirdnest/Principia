#pragma once

#include <concepts>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/tags.hpp"
#include "numerics/concepts.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/matrix_views.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _unbounded_arrays {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::numerics::_concepts;
using namespace principia::numerics::_fixed_arrays;
using namespace principia::numerics::_matrix_views;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

// An allocator that does not initialize the allocated objects.
template<class T>
class uninitialized_allocator : public std::allocator<T> {
 public:
  template<class U, class... Args>
  void construct(U* p, Args&&... args);
};

template<typename Scalar>
class UnboundedMatrix;
template<typename Scalar>
class UnboundedUpperTriangularMatrix;

// The following classes are similar to those in fixed_arrays.hpp, but they have
// an Extend method to add more entries to the arrays.

template<typename Scalar_>
class UnboundedVector final {
 public:
  using Scalar = Scalar_;

  explicit UnboundedVector(std::int64_t size);  // Zero-initialized.
  UnboundedVector(std::int64_t size, uninitialized_t);
  UnboundedVector(std::initializer_list<Scalar> data);

  template<std::int64_t size_>
  explicit UnboundedVector(FixedVector<Scalar, size_> const& data);

  // Constructs an unbounded vector by copying data from the view.  Note that
  // the result is unbounded even if the matrix being viewed is a FixedMatrix.
  template<typename T>
    requires std::same_as<typename T::Scalar, Scalar_>
  explicit UnboundedVector(ColumnView<T> const& view);

  friend bool operator==(UnboundedVector const& left,
                         UnboundedVector const& right) = default;
  friend bool operator!=(UnboundedVector const& left,
                         UnboundedVector const& right) = default;

  Scalar& operator[](std::int64_t index);
  Scalar const& operator[](std::int64_t index) const;

  UnboundedVector& operator=(std::initializer_list<Scalar> right);

  UnboundedVector& operator+=(UnboundedVector const& right);
  UnboundedVector& operator-=(UnboundedVector const& right);
  UnboundedVector& operator*=(double right);
  UnboundedVector& operator/=(double right);

  void Extend(std::int64_t extra_size);
  void Extend(std::int64_t extra_size, uninitialized_t);
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(std::int64_t begin_index);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  UnboundedVector<double> Normalize() const;

  std::int64_t size() const;

  typename std::vector<Scalar>::const_iterator begin() const;
  typename std::vector<Scalar>::const_iterator end() const;

  template<typename H>
  friend H AbslHashValue(H h, UnboundedVector const& vector) {
    for (auto const& scalar : vector) {
      h = H::combine(std::move(h), scalar / si::Unit<Scalar>);
    }
    return h;
  }

 private:
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;
};

template<typename Scalar_>
class UnboundedMatrix final {
 public:
  using Scalar = Scalar_;

  UnboundedMatrix(std::int64_t rows, std::int64_t columns);
  UnboundedMatrix(std::int64_t rows, std::int64_t columns, uninitialized_t);

  // The `data` must be in row-major format and must be for a square matrix.
  UnboundedMatrix(std::initializer_list<Scalar> data);

  UnboundedMatrix(std::int64_t rows, std::int64_t columns,
                  std::initializer_list<Scalar> data);

  explicit UnboundedMatrix(TransposedView<UnboundedMatrix<Scalar>> const& view);

  friend bool operator==(UnboundedMatrix const& left,
                         UnboundedMatrix const& right) = default;
  friend bool operator!=(UnboundedMatrix const& left,
                         UnboundedMatrix const& right) = default;

  // For  0 ≤ i < rows and 0 ≤ j < columns, the entry a_ij is accessed as
  // `a(i, j)`.  If i and j do not satisfy these conditions, the expression
  // `a(i, j)` implies undefined behaviour.
  Scalar& operator()(std::int64_t row, std::int64_t column);
  Scalar const& operator()(std::int64_t row, std::int64_t column) const;

  // Applies the matrix as a bilinear form.  Present for compatibility with
  // `SymmetricBilinearForm`.  Prefer to use `TransposedView` and `operator*`.
  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>>
      operator()(UnboundedVector<LScalar> const& left,
                 UnboundedVector<RScalar> const& right) const;

  UnboundedMatrix& operator=(std::initializer_list<Scalar> right);

  UnboundedMatrix& operator+=(UnboundedMatrix const& right);
  UnboundedMatrix& operator-=(UnboundedMatrix const& right);
  UnboundedMatrix& operator*=(double right);
  UnboundedMatrix& operator/=(double right);

  UnboundedMatrix& operator*=(UnboundedMatrix<double> const& right);

  std::int64_t rows() const;
  std::int64_t columns() const;

  Scalar FrobeniusNorm() const;

  static UnboundedMatrix Identity(std::int64_t rows, std::int64_t columns);

 private:
  std::int64_t rows_;
  std::int64_t columns_;
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(std::ostream& out,
                                  UnboundedMatrix<S> const& matrix);
};

template<typename Scalar_>
class UnboundedLowerTriangularMatrix final {
 public:
  using Scalar = Scalar_;

  explicit UnboundedLowerTriangularMatrix(std::int64_t rows);
  UnboundedLowerTriangularMatrix(std::int64_t rows, uninitialized_t);

  // The `data` must be in row-major format.
  UnboundedLowerTriangularMatrix(std::initializer_list<Scalar> data);

  explicit UnboundedLowerTriangularMatrix(
      TransposedView<UnboundedUpperTriangularMatrix<Scalar>> const& view);

  explicit operator UnboundedMatrix<Scalar>() const;

  friend bool operator==(UnboundedLowerTriangularMatrix const& left,
                         UnboundedLowerTriangularMatrix const& right) = default;
  friend bool operator!=(UnboundedLowerTriangularMatrix const& left,
                         UnboundedLowerTriangularMatrix const& right) = default;

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as `a(i, j)`.
  // If i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  Scalar& operator()(std::int64_t row, std::int64_t column);
  Scalar const& operator()(std::int64_t row, std::int64_t column) const;

  UnboundedLowerTriangularMatrix& operator=(
      std::initializer_list<Scalar> right);

  void Extend(std::int64_t extra_rows);
  void Extend(std::int64_t extra_rows, uninitialized_t);

  // The `data` must be in row-major format.
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(std::int64_t begin_row_index);

  std::int64_t rows() const;
  std::int64_t columns() const;

 private:
  std::int64_t rows_;
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(
      std::ostream& out,
      UnboundedLowerTriangularMatrix<S> const& matrix);
};

template<typename Scalar_>
class UnboundedStrictlyUpperTriangularMatrix final {
 public:
  using Scalar = Scalar_;

  explicit UnboundedStrictlyUpperTriangularMatrix(std::int64_t columns);
  UnboundedStrictlyUpperTriangularMatrix(std::int64_t columns, uninitialized_t);

  // The `data` must be in row-major format.
  UnboundedStrictlyUpperTriangularMatrix(
      std::initializer_list<Scalar> const& data);

  explicit UnboundedStrictlyUpperTriangularMatrix(
      TransposedView<UnboundedLowerTriangularMatrix<Scalar>> const& view);

  explicit operator UnboundedMatrix<Scalar>() const;

  friend bool operator==(
      UnboundedStrictlyUpperTriangularMatrix const& left,
      UnboundedStrictlyUpperTriangularMatrix const& right) = default;
  friend bool operator!=(
      UnboundedStrictlyUpperTriangularMatrix const& left,
      UnboundedStrictlyUpperTriangularMatrix const& right) = default;

  // For  0 ≤ i < j < columns, the entry a_ij is accessed as `a(i, j)`.
  // If i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  Scalar& operator()(std::int64_t row, std::int64_t column);
  Scalar const& operator()(std::int64_t row, std::int64_t column) const;

  UnboundedStrictlyUpperTriangularMatrix& operator=(
      std::initializer_list<Scalar> right);

  void Extend(std::int64_t extra_columns);
  void Extend(std::int64_t extra_columns, uninitialized_t);

  // The `data` must be in row-major format.
  void Extend(std::initializer_list<Scalar> const& data);

  void EraseToEnd(std::int64_t begin_column_index);

  std::int64_t rows() const;
  std::int64_t columns() const;

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This translates a trapezoidal slice to make it column-major.
  static std::vector<Scalar, uninitialized_allocator<Scalar>> Transpose(
      std::initializer_list<Scalar> const& data,
      std::int64_t current_columns,
      std::int64_t extra_columns);

  std::int64_t columns_;
  // Stored in column-major format, so the data passed the public API must be
  // transposed.
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(
      std::ostream& out,
      UnboundedStrictlyUpperTriangularMatrix<S> const& matrix);

  template<typename R>
  friend class Row;
};

template<typename Scalar_>
class UnboundedUpperTriangularMatrix final {
 public:
  using Scalar = Scalar_;

  explicit UnboundedUpperTriangularMatrix(std::int64_t columns);
  UnboundedUpperTriangularMatrix(std::int64_t columns, uninitialized_t);

  // The `data` must be in row-major format.
  UnboundedUpperTriangularMatrix(std::initializer_list<Scalar> const& data);

  explicit UnboundedUpperTriangularMatrix(
      TransposedView<UnboundedLowerTriangularMatrix<Scalar>> const& view);

  explicit operator UnboundedMatrix<Scalar>() const;

  friend bool operator==(UnboundedUpperTriangularMatrix const& left,
                         UnboundedUpperTriangularMatrix const& right) = default;
  friend bool operator!=(UnboundedUpperTriangularMatrix const& left,
                         UnboundedUpperTriangularMatrix const& right) = default;

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as `a(i, j)`.
  // If i and j do not satisfy these conditions, the expression `a(i, j)`
  // implies undefined behaviour.
  Scalar& operator()(std::int64_t row, std::int64_t column);
  Scalar const& operator()(std::int64_t row, std::int64_t column) const;

  UnboundedUpperTriangularMatrix& operator=(
      std::initializer_list<Scalar> right);

  void Extend(std::int64_t extra_columns);
  void Extend(std::int64_t extra_columns, uninitialized_t);

  // The `data` must be in row-major format.
  void Extend(std::initializer_list<Scalar> const& data);

  void EraseToEnd(std::int64_t begin_column_index);

  std::int64_t rows() const;
  std::int64_t columns() const;

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This translates a trapezoidal slice to make it column-major.
  static std::vector<Scalar, uninitialized_allocator<Scalar>> Transpose(
      std::initializer_list<Scalar> const& data,
      std::int64_t current_columns,
      std::int64_t extra_columns);

  std::int64_t columns_;
  // Stored in column-major format, so the data passed the public API must be
  // transposed.
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(
      std::ostream& out,
      UnboundedUpperTriangularMatrix<S> const& matrix);

  template<typename R>
  friend class Row;
};

// Prefer using the operator* that takes a TransposedView.
template<typename LScalar, typename RScalar>
Product<LScalar, RScalar> InnerProduct(
    UnboundedVector<LScalar> const& left,
    UnboundedVector<RScalar> const& right);

template<typename Scalar>
UnboundedVector<double> Normalize(UnboundedVector<Scalar> const& vector);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Product<LScalar, RScalar>> SymmetricProduct(
    UnboundedVector<LScalar> const& left,
    UnboundedVector<RScalar> const& right);

template<typename Scalar>
UnboundedMatrix<Square<Scalar>> SymmetricSquare(
    UnboundedVector<Scalar> const& vector);

// Additive groups.

template<typename Scalar>
UnboundedVector<Scalar> operator+(UnboundedVector<Scalar> const& right);

template<typename Scalar>
UnboundedMatrix<Scalar> operator+(UnboundedMatrix<Scalar> const& right);

template<typename Scalar>
UnboundedVector<Scalar> operator-(UnboundedVector<Scalar> const& right);

template<typename Scalar>
UnboundedMatrix<Scalar> operator-(UnboundedMatrix<Scalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedVector<Sum<LScalar, RScalar>> operator+(
    UnboundedVector<LScalar> const& left,
    UnboundedVector<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Sum<LScalar, RScalar>> operator+(
    UnboundedMatrix<LScalar> const& left,
    UnboundedMatrix<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedVector<Difference<LScalar, RScalar>> operator-(
    UnboundedVector<LScalar> const& left,
    UnboundedVector<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Difference<LScalar, RScalar>> operator-(
    UnboundedMatrix<LScalar> const& left,
    UnboundedMatrix<RScalar> const& right);

// Vector spaces.

template<typename LScalar, typename RScalar>
UnboundedVector<Product<LScalar, RScalar>> operator*(
    LScalar const& left,
    UnboundedVector<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedVector<Product<LScalar, RScalar>> operator*(
    UnboundedVector<LScalar> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Product<LScalar, RScalar>>
operator*(LScalar const& left,
          UnboundedMatrix<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Product<LScalar, RScalar>>
operator*(UnboundedMatrix<LScalar> const& left,
          RScalar const& right);

template<typename LScalar, typename RScalar>
UnboundedVector<Quotient<LScalar, RScalar>> operator/(
    UnboundedVector<LScalar> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Quotient<LScalar, RScalar>>
operator/(UnboundedMatrix<LScalar> const& left,
          RScalar const& right);

// Hilbert space and algebra.

// TODO(phl): fixed_arrays.hpp has an operator* that takes a row.

template<typename LScalar, typename RScalar>
Product<LScalar, RScalar> operator*(
    TransposedView<UnboundedVector<LScalar>> const& left,
    UnboundedVector<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Product<LScalar, RScalar>> operator*(
    UnboundedVector<LScalar> const& left,
    TransposedView<UnboundedVector<RScalar>> const& right);

template<typename LScalar, typename RScalar>
UnboundedMatrix<Product<LScalar, RScalar>> operator*(
    UnboundedMatrix<LScalar> const& left,
    UnboundedMatrix<RScalar> const& right);

template<typename LScalar, typename RScalar>
UnboundedVector<Product<LScalar, RScalar>> operator*(
    UnboundedMatrix<LScalar> const& left,
    UnboundedVector<RScalar> const& right);

template<typename LMatrix, typename RScalar>
UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> operator*(
    BlockView<LMatrix> const& left,
    UnboundedVector<RScalar> const& right);

template<typename LMatrix, typename RScalar>
UnboundedVector<Product<typename LMatrix::Scalar, RScalar>> operator*(
    TransposedView<BlockView<LMatrix>> const& left,
    UnboundedVector<RScalar> const& right);

// Use this operator to multiply a row vector with a matrix.  We don't have an
// operator returning a TransposedView as that would cause dangling references.
template<typename LScalar, typename RScalar>
UnboundedVector<Product<LScalar, RScalar>> operator*(
    TransposedView<UnboundedMatrix<LScalar>> const& left,
    UnboundedVector<RScalar> const& right);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedVector<Scalar> const& vector);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedMatrix<Scalar> const& matrix);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedLowerTriangularMatrix<Scalar> const& matrix);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedUpperTriangularMatrix<Scalar> const& matrix);

}  // namespace internal

using internal::UnboundedLowerTriangularMatrix;
using internal::UnboundedMatrix;
using internal::UnboundedStrictlyUpperTriangularMatrix;
using internal::UnboundedUpperTriangularMatrix;
using internal::UnboundedVector;

}  // namespace _unbounded_arrays
}  // namespace numerics
}  // namespace principia

#include "numerics/unbounded_arrays_body.hpp"
