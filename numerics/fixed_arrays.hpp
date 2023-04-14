#pragma once

#include <array>
#include <vector>

#include "base/tags.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _fixed_arrays {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::numerics::_transposed_view;
using namespace principia::quantities::_named_quantities;

template<typename Scalar, int rows, int columns>
class FixedMatrix;

template<typename Scalar, int size_>
class FixedVector final {
 public:
  constexpr FixedVector();
  explicit FixedVector(uninitialized_t);

  // TODO(egg): Figure out why we have a move-conversion for |FixedVector| but
  // not the matrices.
  constexpr FixedVector(std::array<Scalar, size_> const& data);
  constexpr FixedVector(std::array<Scalar, size_>&& data);

  TransposedView<FixedVector> Transpose() const;

  Scalar Norm() const;

  static constexpr int size() { return size_; }

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  explicit operator std::vector<Scalar>() const;

  bool operator==(FixedVector const& right) const;
  bool operator!=(FixedVector const& right) const;

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
  constexpr FixedMatrix(std::array<Scalar, size()> const& data);

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

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr FixedVector<Quotient<ScalarLeft, ScalarRight>, size> operator/(
    FixedVector<ScalarLeft, size> const& left,
    ScalarRight const& right);

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr FixedVector<Difference<ScalarLeft, ScalarRight>, size> operator-(
    FixedVector<ScalarLeft, size> const& left,
    FixedVector<ScalarRight, size> const& right);

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr Product<ScalarLeft, ScalarRight> operator*(
    ScalarLeft* const left,
    FixedVector<ScalarRight, size> const& right);

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr Product<ScalarLeft, ScalarRight> operator*(
    TransposedView<FixedVector<ScalarLeft, size>> const& left,
    FixedVector<ScalarRight, size> const& right);

template<typename ScalarLeft, typename ScalarRight,
         int rows, int dimension, int columns>
constexpr FixedMatrix<Product<ScalarLeft, ScalarRight>, rows, columns>
operator*(FixedMatrix<ScalarLeft, rows, dimension> const& left,
          FixedMatrix<ScalarRight, dimension, columns> const& right);

template<typename ScalarLeft, typename ScalarRight, int rows, int columns>
constexpr FixedVector<Product<ScalarLeft, ScalarRight>, rows> operator*(
    FixedMatrix<ScalarLeft, rows, columns> const& left,
    FixedVector<ScalarRight, columns> const& right);

template<typename Scalar, int size>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size> const& vector);

template<typename Scalar, int rows, int columns>
std::ostream& operator<<(std::ostream& out,
                         FixedMatrix<Scalar, rows, columns> const& matrix);

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
