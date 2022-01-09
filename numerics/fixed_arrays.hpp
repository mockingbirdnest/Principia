
#pragma once

#include <array>
#include <vector>

#include "base/tags.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_fixed_arrays {

using base::uninitialized_t;
using quantities::Difference;
using quantities::Product;

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

  bool operator==(FixedVector const& right) const;

  constexpr Scalar& operator[](int index);
  constexpr Scalar const& operator[](int index) const;

  explicit operator std::vector<Scalar>() const;

  static constexpr int size = size_;

 private:
  std::array<Scalar, size> data_;

  template<typename L, typename R, int r, int c>
  friend FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
  template<typename S, int s>
  friend std::ostream& operator<<(std::ostream& out,
                                  FixedVector<S, s> const& vector);
};

template<typename Scalar, int rows, int columns>
class FixedMatrix final {
 public:
  constexpr FixedMatrix();
  explicit FixedMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedMatrix(std::array<Scalar, rows * columns> const& data);

  bool operator==(FixedMatrix const& right) const;

  template<int r>
  class Row {
   public:
    explicit Row(const FixedMatrix* matrix);

    constexpr Scalar const& operator[](int index) const;

    // The template deduction runs into trouble if this operator is declared at
    // namespace scope.
    template<typename S>
    Product<Scalar, S> operator*(FixedVector<S, columns> const& right);

   private:
    const FixedMatrix* matrix_;
  };

  template<int r>
  typename FixedMatrix::template Row<r> row() const;

  // For  0 < i < rows and 0 < j < columns, the entry a_ij is accessed as
  // |a[i][j]|.  if i and j do not satisfy these conditions, the expression
  // |a[i][j]| is erroneous.
  Scalar* operator[](int index);
  constexpr Scalar const* operator[](int index) const;

 private:
  std::array<Scalar, rows * columns> data_;

  template<typename L, typename R, int r, int c>
  friend FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
};

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr FixedVector<Difference<ScalarLeft, ScalarRight>, size> operator-(
    FixedVector<ScalarLeft, size> const& left,
    FixedVector<ScalarRight, size> const& right);

template<typename ScalarLeft, typename ScalarRight, int rows, int columns>
FixedVector<Product<ScalarLeft, ScalarRight>, rows> operator*(
    FixedMatrix<ScalarLeft, rows, columns> const& left,
    FixedVector<ScalarRight, columns> const& right);

template<typename Scalar, int rows_>
class FixedStrictlyLowerTriangularMatrix final {
 public:
  static constexpr int size = rows_;
  static constexpr int dimension = rows_ * (rows_ - 1) / 2;

  constexpr FixedStrictlyLowerTriangularMatrix();
  explicit FixedStrictlyLowerTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedStrictlyLowerTriangularMatrix(
      std::array<Scalar, dimension> const& data);

  bool operator==(FixedStrictlyLowerTriangularMatrix const& right) const;

  // For  0 < j < i < rows, the entry a_ij is accessed as |a[i][j]|.
  // if i and j do not satisfy these conditions, the expression |a[i][j]| is
  // erroneous.
  Scalar* operator[](int index);
  constexpr Scalar const* operator[](int index) const;

 private:
  std::array<Scalar, dimension> data_;
};

template<typename Scalar, int rows_>
class FixedLowerTriangularMatrix final {
 public:
  static constexpr int rows = rows_;
  static constexpr int dimension = rows * (rows + 1) / 2;

  constexpr FixedLowerTriangularMatrix();
  explicit FixedLowerTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedLowerTriangularMatrix(
      std::array<Scalar, dimension> const& data);

  bool operator==(FixedLowerTriangularMatrix const& right) const;

  // For  0 < j <= i < rows, the entry a_ij is accessed as |a[i][j]|.
  // if i and j do not satisfy these conditions, the expression |a[i][j]| is
  // erroneous.
  Scalar* operator[](int index);
  constexpr Scalar const* operator[](int index) const;

 private:
  std::array<Scalar, dimension> data_;
};

template<typename Scalar, int columns_>
class FixedUpperTriangularMatrix final {
 public:
  static constexpr int columns() { return columns_; }
  static constexpr int dimension = columns_ * (columns_ + 1) / 2;

  constexpr FixedUpperTriangularMatrix();
  explicit FixedUpperTriangularMatrix(uninitialized_t);

  // The |data| must be in row-major format.
  constexpr FixedUpperTriangularMatrix(
      std::array<Scalar, dimension> const& data);

  bool operator==(FixedUpperTriangularMatrix const& right) const;

  // A helper class for indexing column-major data in a human-friendly manner.
  template<typename Matrix>
  class Row {
   public:
    Scalar& operator[](int column);
    Scalar const& operator[](int column) const;

   private:
    explicit Row(Matrix& matrix, int row);

    // We need to remove the const because, when this class is instantiated with
    // 'FixedUpperTriangularMatrix const', the first operator[], not the second,
    // is picked by overload resolution.
    std::remove_const_t<Matrix>& matrix_;
    int row_;

    template<typename S, int c>
    friend class FixedUpperTriangularMatrix;
  };

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as |a[i][j]|.
  // if i and j do not satisfy these conditions, the expression |a[i][j]| is
  // erroneous.
  Row<FixedUpperTriangularMatrix> operator[](int row);
  Row<FixedUpperTriangularMatrix const> operator[](int row) const;

 private:
  // For ease of writing matrices in tests, the input data is received in row-
  // major format.  This translates a trapezoidal slice to make it column-major.
  static std::array<Scalar, dimension> Transpose(
      std::array<Scalar, dimension> const& data);

  std::array<Scalar, dimension> data_;
};

template<typename Scalar, int size>
std::ostream& operator<<(std::ostream& out,
                         FixedVector<Scalar, size> const& vector);

template<typename Scalar, int rows>
std::ostream& operator<<(
    std::ostream& out,
    FixedLowerTriangularMatrix<Scalar, rows> const& matrix);

template<typename Scalar, int columns>
std::ostream& operator<<(
    std::ostream& out,
    FixedUpperTriangularMatrix<Scalar, columns> const& matrix);

}  // namespace internal_fixed_arrays

using internal_fixed_arrays::FixedLowerTriangularMatrix;
using internal_fixed_arrays::FixedMatrix;
using internal_fixed_arrays::FixedStrictlyLowerTriangularMatrix;
using internal_fixed_arrays::FixedUpperTriangularMatrix;
using internal_fixed_arrays::FixedVector;

}  // namespace numerics
}  // namespace principia

#include "numerics/fixed_arrays_body.hpp"
