
#pragma once

#include <initializer_list>
#include <memory>
#include <type_traits>
#include <vector>

#include "base/tags.hpp"

namespace principia {
namespace numerics {
namespace internal_unbounded_arrays {

using base::uninitialized_t;

// An allocator that does not initialize the allocated objects.
template<class T>
class uninitialized_allocator : public std::allocator<T> {
 public:
  template<class U, class... Args>
  void construct(U* p, Args&&... args);
};

template<typename Scalar>
class UnboundedUpperTriangularMatrix;

// The following classes are similar to those in fixed_arrays.hpp, but they have
// an Extend method to add more entries to the arrays.

template<typename Scalar>
class UnboundedVector final {
 public:
  explicit UnboundedVector(int size);
  UnboundedVector(int size, uninitialized_t);
  UnboundedVector(std::initializer_list<Scalar> data);

  void Extend(int extra_size);
  void Extend(int extra_size, uninitialized_t);
  void Extend(std::initializer_list<Scalar> data);

  void EraseToEnd(int begin_index);

  int size() const;

  bool operator==(UnboundedVector const& right) const;

  Scalar& operator[](int index);
  Scalar const& operator[](int index) const;

 private:
  std::vector<Scalar, uninitialized_allocator<Scalar>> data_;

  template<typename S>
  friend std::ostream& operator<<(std::ostream& out,
                                  UnboundedVector<S> const& vector);
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

  UnboundedUpperTriangularMatrix<Scalar> Transpose() const;

  int rows() const;
  int dimension() const;

  bool operator==(UnboundedLowerTriangularMatrix const& right) const;

  // For  0 ≤ j ≤ i < rows, the entry a_ij is accessed as |a[i][j]|.
  // if i and j do not satisfy these conditions, the expression |a[i][j]| is
  // erroneous.
  Scalar* operator[](int index);
  Scalar const* operator[](int index) const;

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

  UnboundedLowerTriangularMatrix<Scalar> Transpose() const;

  int columns() const;
  int dimension() const;

  bool operator==(UnboundedUpperTriangularMatrix const& right) const;

  // A helper class for indexing column-major data in a human-friendly manner.
  template<typename Matrix>
  class Row {
   public:
    Scalar& operator[](int column);
    Scalar const& operator[](int column) const;

   private:
    explicit Row(Matrix& matrix, int row);

    // We need to remove the const because, when this class is instantiated with
    // 'UnboundedUpperTriangularMatrix const', the first operator[], not the
    // second, is picked by overload resolution.
    std::remove_const_t<Matrix>& matrix_;
    int row_;

    template<typename S>
    friend class UnboundedUpperTriangularMatrix;
  };

  // For  0 ≤ i ≤ j < columns, the entry a_ij is accessed as |a[i][j]|.
  // if i and j do not satisfy these conditions, the expression |a[i][j]| is
  // erroneous.
  Row<UnboundedUpperTriangularMatrix> operator[](int row);
  Row<UnboundedUpperTriangularMatrix const> operator[](int row) const;

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

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedVector<Scalar> const& vector);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedLowerTriangularMatrix<Scalar> const& matrix);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         UnboundedUpperTriangularMatrix<Scalar> const& matrix);

}  // namespace internal_unbounded_arrays

using internal_unbounded_arrays::UnboundedLowerTriangularMatrix;
using internal_unbounded_arrays::UnboundedUpperTriangularMatrix;
using internal_unbounded_arrays::UnboundedVector;

}  // namespace numerics
}  // namespace principia

#include "numerics/unbounded_arrays_body.hpp"
