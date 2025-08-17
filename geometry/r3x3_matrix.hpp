#pragma once

#include <array>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

#include "base/not_null.hpp"
#include "base/tags.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/arithmetic.hpp"
#include "quantities/concepts.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _r3x3_matrix {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_tags;
using namespace principia::geometry::_r3_element;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_concepts;

// An `R3x3Matrix` is an element of the associative algebra of 3-by-3 matrices
// over `Scalar`.  `Scalar` should be a vector space over ℝ, represented by
// `double`.
template<typename Scalar>
class R3x3Matrix final {
 public:
  R3x3Matrix() = default;
  explicit R3x3Matrix(uninitialized_t);
  R3x3Matrix(R3Element<Scalar> const& row_x,
             R3Element<Scalar> const& row_y,
             R3Element<Scalar> const& row_z);

  friend bool operator==(R3x3Matrix const& left,
                         R3x3Matrix const& right) = default;
  friend bool operator!=(R3x3Matrix const& left,
                         R3x3Matrix const& right) = default;

  Scalar operator()(int r, int c) const;
  Scalar& operator()(int r, int c);

  R3x3Matrix& operator+=(R3x3Matrix const& right);
  R3x3Matrix& operator-=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(double right);
  R3x3Matrix& operator/=(double right);

  static R3x3Matrix DiagonalMatrix(R3Element<Scalar> const& diagonal);

  R3Element<Scalar> Diagonal() const;

  Scalar Trace() const;
  Cube<Scalar> Determinant() const;
  R3x3Matrix Transpose() const;

  Scalar FrobeniusNorm() const;

  void QRDecomposition(R3x3Matrix<double>& Q, R3x3Matrix& R) const;

  template<typename RScalar>
  R3Element<Quotient<RScalar, Scalar>> Solve(
      R3Element<RScalar> const& rhs) const;

  R3Element<Scalar> const& row_x() const;
  R3Element<Scalar> const& row_y() const;
  R3Element<Scalar> const& row_z() const;

  template<typename S = Scalar,
           typename = std::enable_if_t<std::is_same_v<S, double>>>
  static R3x3Matrix<S> Identity();

  void WriteToMessage(not_null<serialization::R3x3Matrix*> message) const;
  static R3x3Matrix ReadFromMessage(serialization::R3x3Matrix const& message);

  // Clang on MacOS don't like all these friends.
#if !OS_MACOSX

 private:
#endif
  enum Indices { X = 0, Y = 1, Z = 2 };
  std::array<R3Element<Scalar>, 3> rows_;

#if !OS_MACOSX
  template<typename S>
  friend class R3x3Matrix;

  template<typename S>
  friend R3x3Matrix<S> operator+(R3x3Matrix<S> const& right);
  template<typename S>
  friend R3x3Matrix<S> operator-(R3x3Matrix<S> const& right);

  template<typename S>
  friend R3x3Matrix<S> operator+(R3x3Matrix<S> const& left,
                                 R3x3Matrix<S> const& right);
  template<typename S>
  friend R3x3Matrix<S> operator-(R3x3Matrix<S> const& left,
                                 R3x3Matrix<S> const& right);
  template<typename LS, typename RS>
  friend R3x3Matrix<Product<LS, RS>> operator*(R3x3Matrix<LS> const& left,
                                               R3x3Matrix<RS> const& right);

  template<typename LS, typename RS>
  friend R3Element<Product<LS, RS>> operator*(R3x3Matrix<LS> const& left,
                                              R3Element<RS> const& right);
  template<typename LS, typename RS>
  friend R3Element<Product<LS, RS>> operator*(R3Element<LS> const& left,
                                              R3x3Matrix<RS> const& right);

  template<typename LS, typename RS>
    requires convertible_to_quantity<LS>
  friend R3x3Matrix<Product<LS, RS>> operator*(LS const& left,
                                               R3x3Matrix<RS> const& right);
  template<typename LS, typename RS>
    requires convertible_to_quantity<RS>
  friend R3x3Matrix<Product<LS, RS>> operator*(R3x3Matrix<LS> const& left,
                                               RS const& right);
  template<typename LS, typename RS>
    requires convertible_to_quantity<RS>
  friend R3x3Matrix<Quotient<LS, RS>> operator/(R3x3Matrix<LS> const& left,
                                                RS const& right);

  template<typename S>
  friend std::string DebugString(R3x3Matrix<S> const& r3x3_matrix);
#endif
};

template<typename Scalar>
R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& right);
template<typename Scalar>
R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& right);

template<typename Scalar>
R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& left,
                             R3x3Matrix<Scalar> const& right);
template<typename Scalar>
R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& left,
                             R3x3Matrix<Scalar> const& right);

template<typename LScalar, typename RScalar>
R3x3Matrix<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    R3x3Matrix<RScalar> const& right);

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    R3Element<RScalar> const& right);
template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    R3Element<LScalar> const& left,
    R3x3Matrix<RScalar> const& right);

template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<LScalar>
R3x3Matrix<Product<LScalar, RScalar>> operator*(
    LScalar const& left,
    R3x3Matrix<RScalar> const& right);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3x3Matrix<Product<LScalar, RScalar>> operator*(R3x3Matrix<LScalar> const& left,
                                                RScalar const& right);
template<typename LScalar, typename RScalar>
  requires convertible_to_quantity<RScalar>
R3x3Matrix<Quotient<LScalar, RScalar>> operator/(
    R3x3Matrix<LScalar> const& left,
    RScalar const& right);

template<typename LScalar, typename RScalar>
R3x3Matrix<Product<LScalar, RScalar>> KroneckerProduct(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right);

template<typename Scalar>
bool operator==(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right);

template<typename Scalar>
std::string DebugString(R3x3Matrix<Scalar> const& r3x3_matrix);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix<Scalar> const& r3x3_matrix);

}  // namespace internal

using internal::KroneckerProduct;
using internal::R3x3Matrix;

}  // namespace _r3x3_matrix
}  // namespace geometry
}  // namespace principia

#include "geometry/r3x3_matrix_body.hpp"
