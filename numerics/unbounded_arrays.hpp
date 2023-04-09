#pragma once

#include <initializer_list>
#include <memory>
#include <type_traits>
#include <vector>

#include "base/tags.hpp"
#include "numerics/transposed_view.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace numerics {
namespace _unbounded_arrays {
namespace internal {

using namespace principia::base::_tags;
using namespace principia::quantities::_named_quantities;

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

template<typename Scalar>
class UnboundedVector final {
 public:
  explicit UnboundedVector(int size);  // Zero-initialized.
  UnboundedVector(int size, uninitialized_t);
  UnboundedVector(std::initializer_list<Scalar> data);

  TransposedView<UnboundedVector> Transpose() const;

  void Extend(int extra_size);
  void Extend(int extra_size, uninitialized_t);
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(int begin_index);

  Scalar Norm() const;

  int size() const;

  Scalar& operator[](int index);
  Scalar const& operator[](int index) const;

  bool operator==(UnboundedVector const& right) const;
  bool operator!=(UnboundedVector const& right) const;

 private:
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;
};

template<typename Scalar>
class UnboundedMatrix final {
 public:
  UnboundedMatrix(int rows, int columns);
  UnboundedMatrix(int rows, int columns, uninitialized_t);

  // The |data| must be in row-major format.
  UnboundedMatrix(std::initializer_list<Scalar> data);

  int columns() const;
  int rows() const;
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

  static UnboundedMatrix Identity(int rows, int columns);

 private:
  int rows_;
  int columns_;
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(std::ostream& out,
                                  UnboundedMatrix<S> const& matrix);
};

template<typename Scalar>
class UnboundedLowerTriangularMatrix final {
 public:
  explicit UnboundedLowerTriangularMatrix(int rows);
  UnboundedLowerTriangularMatrix(int rows, uninitialized_t);

  // The |data| must be in row-major format.
  UnboundedLowerTriangularMatrix(std::initializer_list<Scalar> data);

  void Extend(int extra_rows);
  void Extend(int extra_rows, uninitialized_t);

  // The |data| must be in row-major format.
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(int begin_row_index);

  int columns() const;
  int rows() const;
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

template<typename Scalar>
class UnboundedUpperTriangularMatrix final {
 public:
  explicit UnboundedUpperTriangularMatrix(int columns);
  UnboundedUpperTriangularMatrix(int columns, uninitialized_t);

  // The |data| must be in row-major format.
  UnboundedUpperTriangularMatrix(std::initializer_list<Scalar> const& data);

  void Extend(int extra_columns);
  void Extend(int extra_columns, uninitialized_t);

  // The |data| must be in row-major format.
  void Extend(std::initializer_list<Scalar> const& data);

  void EraseToEnd(int begin_column_index);

  int columns() const;
  int rows() const;
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

template<typename ScalarLeft, typename ScalarRight>
UnboundedVector<Quotient<ScalarLeft, ScalarRight>> operator/(
    UnboundedVector<ScalarLeft> const& left,
    ScalarRight const& right);

template<typename ScalarLeft, typename ScalarRight>
Product<ScalarLeft, ScalarRight> operator*(
    TransposedView<UnboundedVector<ScalarLeft>> const& left,
    UnboundedVector<ScalarRight> const& right);

template<typename ScalarLeft, typename ScalarRight>
UnboundedMatrix<Product<ScalarLeft, ScalarRight>> operator*(
    UnboundedMatrix<ScalarLeft> const& left,
    UnboundedMatrix<ScalarRight> const& right);

template<typename ScalarLeft, typename ScalarRight>
UnboundedVector<Product<ScalarLeft, ScalarRight>> operator*(
    UnboundedMatrix<ScalarLeft> const& left,
    UnboundedVector<ScalarRight> const& right);

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

namespace principia::numerics {
using namespace principia::numerics::_unbounded_arrays;
}  // namespace principia::numerics

#include "numerics/unbounded_arrays_body.hpp"
