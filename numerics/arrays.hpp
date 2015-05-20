#pragma once

#include <array>
#include <initializer_list>

namespace principia {
namespace numerics {

template<typename Scalar, int rows, int columns>
class FixedMatrix;

template<typename Scalar, int size>
class FixedVector {
 public:
  FixedVector();  // Zero-initialized.
  FixedVector(std::array<Scalar, size> const& data);
  FixedVector(std::initializer_list<Scalar> const& data);

  bool operator==(FixedVector const& right) const;

  FixedVector& operator=(std::initializer_list<Scalar> const& right);

 private:
  std::array<Scalar, size> data_;

  template<typename S, int r, int c>
  friend FixedVector<S, r> operator*(
      FixedMatrix<S, r, c> const& left,
      FixedVector<S, c> const& right);
};

template<typename Scalar, int rows, int columns>
class FixedMatrix {
 public:
  // The |data| must be in row-major format.
  FixedMatrix(std::array<Scalar, rows * columns> const& data);
  FixedMatrix(std::initializer_list<Scalar> const& data);

  bool operator==(FixedMatrix const& right) const;
  FixedMatrix& operator=(std::initializer_list<Scalar> const& right);

 private:
  std::array<Scalar, rows * columns> data_;

  template<typename Scalar, int r, int c>
  friend FixedVector<Scalar, r> operator*(
      FixedMatrix<Scalar, r, c> const& left,
      FixedVector<Scalar, c> const& right);
};

template<typename Scalar, int rows, int columns>
FixedVector<Scalar, rows> operator*(
    FixedMatrix<Scalar, rows, columns> const& left,
    FixedVector<Scalar, columns> const& right);

}  // namespace numerics
}  // namespace principia

#include "numerics/arrays_body.hpp"
