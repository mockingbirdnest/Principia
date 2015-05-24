#pragma once

#include <array>
#include <initializer_list>
#include <vector>

#include "quantities/quantities.hpp"

namespace principia {

using quantities::Product;

namespace numerics {

template<typename Scalar, int rows, int columns>
class FixedMatrix;

template<typename Scalar, int size>
class FixedVector {
 public:
  FixedVector();  // Zero-initialized.
  explicit FixedVector(std::array<Scalar, size> const& data);
  FixedVector(
      std::initializer_list<Scalar> const& data);  // NOLINT(runtime/explicit)

  bool operator==(FixedVector const& right) const;
  FixedVector& operator=(std::initializer_list<Scalar> const& right);

  Scalar& operator[](int const index);
  Scalar const& operator[](int const index) const;

  operator std::vector<Scalar>() const;

 private:
  std::array<Scalar, size> data_;

  template<typename L, typename R, int r, int c>
  friend FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
};

template<typename Scalar, int rows, int columns>
class FixedMatrix {
 public:
  // The |data| must be in row-major format.
  explicit FixedMatrix(std::array<Scalar, rows * columns> const& data);
  FixedMatrix(
      std::initializer_list<Scalar> const& data);  // NOLINT(runtime/explicit)

  bool operator==(FixedMatrix const& right) const;
  FixedMatrix& operator=(std::initializer_list<Scalar> const& right);

 private:
  std::array<Scalar, rows * columns> data_;

  template<typename L, typename R, int r, int c>
  friend FixedVector<Product<L, R>, r> operator*(
      FixedMatrix<L, r, c> const& left,
      FixedVector<R, c> const& right);
};

template<typename ScalarLeft, typename ScalarRight, int rows, int columns>
FixedVector<Product<ScalarLeft, ScalarRight>, rows> operator*(
    FixedMatrix<ScalarLeft, rows, columns> const& left,
    FixedVector<ScalarRight, columns> const& right);

}  // namespace numerics
}  // namespace principia

#include "numerics/arrays_body.hpp"
