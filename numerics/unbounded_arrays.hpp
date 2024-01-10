#pragma once

#include <initializer_list>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/tags.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _unbounded_arrays {
namespace internal {

// TODO(phl): This should support the same operations as fixed_arrays.hpp.

using namespace principia::base::_tags;
using namespace principia::numerics::_fixed_arrays;
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

  explicit UnboundedVector(int size);  // Zero-initialized.
  UnboundedVector(int size, uninitialized_t);
  UnboundedVector(std::initializer_list<Scalar> data);
  template<int size_>
  explicit UnboundedVector(FixedVector<Scalar, size_> const& data);

  void Extend(int extra_size);
  void Extend(int extra_size, uninitialized_t);
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(int begin_index);

  Scalar Norm() const;
  Square<Scalar> Norm²() const;

  int size() const;

  typename std::vector<Scalar>::const_iterator begin() const;
  typename std::vector<Scalar>::const_iterator end() const;

  Scalar& operator[](int index);
  Scalar const& operator[](int index) const;

  bool operator==(UnboundedVector const& right) const;
  bool operator!=(UnboundedVector const& right) const;

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

  UnboundedMatrix(int rows, int columns);
  UnboundedMatrix(int rows, int columns, uninitialized_t);

  // The |data| must be in row-major format and must be for a square matrix.
  UnboundedMatrix(std::initializer_list<Scalar> data);

  UnboundedMatrix(int rows, int columns, std::initializer_list<Scalar> data);

  int rows() const;
  int columns() const;
  int size() const;

  // For  0 ≤ i < rows and 0 ≤ j < columns, the entry a_ij is accessed as
  // |a(i, j)|.  If i and j do not satisfy these conditions, the expression
  // |a(i, j)| implies undefined behaviour.
  Scalar& operator()(int row, int column);
  Scalar const& operator()(int row, int column) const;

  UnboundedMatrix Transpose() const;

  Scalar FrobeniusNorm() const;

  bool operator==(UnboundedMatrix const& right) const;
  bool operator!=(UnboundedMatrix const& right) const;

  // Applies the matrix as a bilinear form.  Present for compatibility with
  // |SymmetricBilinearForm|.  Prefer to use |TransposedView| and |operator*|.
  template<typename LScalar, typename RScalar>
  Product<Scalar, Product<LScalar, RScalar>>
      operator()(UnboundedVector<LScalar> const& left,
                 UnboundedVector<RScalar> const& right) const;

  static UnboundedMatrix Identity(int rows, int columns);

 private:
  int rows_;
  int columns_;
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(std::ostream& out,
                                  UnboundedMatrix<S> const& matrix);
};

template<typename Scalar_>
class UnboundedLowerTriangularMatrix final {
 public:
  using Scalar = Scalar_;

  explicit UnboundedLowerTriangularMatrix(int rows);
  UnboundedLowerTriangularMatrix(int rows, uninitialized_t);

  // The |data| must be in row-major format.
  UnboundedLowerTriangularMatrix(std::initializer_list<Scalar> data);

  void Extend(int extra_rows);
  void Extend(int extra_rows, uninitialized_t);

  // The |data| must be in row-major format.
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(int begin_row_index);

  int rows() const;
  int columns() const;
  int size() const;

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as |a(i, j)|.
  // If i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  Scalar& operator()(int row, int column);
  Scalar const& operator()(int row, int column) const;

  UnboundedUpperTriangularMatrix<Scalar> Transpose() const;

  bool operator==(UnboundedLowerTriangularMatrix const& right) const;
  bool operator!=(UnboundedLowerTriangularMatrix const& right) const;

 private:
  int rows_;
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(
      std::ostream& out,
      UnboundedLowerTriangularMatrix<S> const& matrix);
};

template<typename Scalar_>
class UnboundedUpperTriangularMatrix final {
 public:
  using Scalar = Scalar_;

  explicit UnboundedUpperTriangularMatrix(int columns);
  UnboundedUpperTriangularMatrix(int columns, uninitialized_t);

  // The |data| must be in row-major format.
  UnboundedUpperTriangularMatrix(std::initializer_list<Scalar> const& data);

  void Extend(int extra_columns);
  void Extend(int extra_columns, uninitialized_t);

  // The |data| must be in row-major format.
  void Extend(std::initializer_list<Scalar> const& data);

  void EraseToEnd(int begin_column_index);

  int rows() const;
  int columns() const;
  int size() const;

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as |a(i, j)|.
  // If i and j do not satisfy these conditions, the expression |a(i, j)|
  // implies undefined behaviour.
  Scalar& operator()(int row, int column);
  Scalar const& operator()(int row, int column) const;

  UnboundedLowerTriangularMatrix<Scalar> Transpose() const;

  bool operator==(UnboundedUpperTriangularMatrix const& right) const;
  bool operator!=(UnboundedUpperTriangularMatrix const& right) const;

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This translates a trapezoidal slice to make it column-major.
  static std::vector<Scalar, uninitialized_allocator<Scalar>> Transpose(
      std::initializer_list<Scalar> const& data,
      int current_columns,
      int extra_columns);

  int columns_;
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
UnboundedVector<Scalar> operator-(
    UnboundedVector<Scalar> const& right);

template<typename Scalar>
UnboundedMatrix<Scalar> operator-(
    UnboundedMatrix<Scalar> const& right);

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

template<typename Scalar>
UnboundedVector<Scalar>& operator+=(
    UnboundedVector<Scalar>& left,
    UnboundedVector<Scalar> const& right);

template<typename Scalar>
UnboundedMatrix<Scalar>& operator+=(
    UnboundedMatrix<Scalar>& left,
    UnboundedMatrix<Scalar> const& right);

template<typename Scalar>
UnboundedVector<Scalar>& operator-=(
    UnboundedVector<Scalar>& left,
    UnboundedVector<Scalar> const& right);

template<typename Scalar>
UnboundedMatrix<Scalar>& operator-=(
    UnboundedMatrix<Scalar>& left,
    UnboundedMatrix<Scalar> const& right);

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

template<typename Scalar>
UnboundedVector<Scalar>& operator*=(
    UnboundedVector<Scalar>& left,
    double right);

template<typename Scalar>
UnboundedMatrix<Scalar>& operator*=(
    UnboundedMatrix<Scalar>& left,
    double right);

template<typename Scalar>
UnboundedVector<Scalar>& operator/=(
    UnboundedVector<Scalar>& left,
    double right);

template<typename Scalar>
UnboundedMatrix<Scalar>& operator/=(
    UnboundedMatrix<Scalar>& left,
    double right);

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
using internal::UnboundedUpperTriangularMatrix;
using internal::UnboundedVector;

}  // namespace _unbounded_arrays
}  // namespace numerics
}  // namespace principia

#include "numerics/unbounded_arrays_body.hpp"
