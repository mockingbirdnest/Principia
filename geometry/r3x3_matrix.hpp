
#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>

#include "geometry/r3_element.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_r3x3_matrix {

using base::not_null;
using quantities::Product;

// An |R3x3Matrix| is an element of the associative algebra of 3-by-3 matrices
// over |Scalar|.  |Scalar| should be a vector space over ℝ, represented by
// |double|.
template<typename Scalar>
class R3x3Matrix final {
 public:
  R3x3Matrix(R3Element<Scalar> const& row_x,
             R3Element<Scalar> const& row_y,
             R3Element<Scalar> const& row_z);

  Scalar Trace() const;
  R3x3Matrix Transpose() const;

  Scalar operator()(int r, int c) const;

  R3x3Matrix& operator+=(R3x3Matrix const& right);
  R3x3Matrix& operator-=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(double right);
  R3x3Matrix& operator/=(double right);

  template<typename = std::is_arithmetic<Scalar>::value>
  static R3x3Matrix<double> Identity();

  void WriteToMessage(not_null<serialization::R3x3Matrix*> message) const;
  static R3x3Matrix ReadFromMessage(serialization::R3x3Matrix const& message);

 private:
  R3Element<Scalar> row_x_;
  R3Element<Scalar> row_y_;
  R3Element<Scalar> row_z_;

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

  template<typename S>
  friend R3x3Matrix<S> operator*(double left, R3x3Matrix<S> const& right);
  template<typename S>
  friend R3x3Matrix<S> operator*(R3x3Matrix<S> const& left, double right);
  template<typename S>
  friend R3x3Matrix<S> operator/(R3x3Matrix<S> const& left, double right);

  template<typename S>
  friend bool operator==(R3x3Matrix<S> const& left,
                         R3x3Matrix<S> const& right);
  template<typename S>
  friend bool operator!=(R3x3Matrix<S> const& left,
                         R3x3Matrix<S> const& right);

  template<typename S>
  friend std::string DebugString(R3x3Matrix<S> const& r3x3_matrix);
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

template<typename Scalar>
R3x3Matrix<Scalar> operator*(double left, R3x3Matrix<Scalar> const& right);
template<typename Scalar>
R3x3Matrix<Scalar> operator*(R3x3Matrix<Scalar> const& left, double right);
template<typename Scalar>
R3x3Matrix<Scalar> operator/(R3x3Matrix<Scalar> const& left, double right);

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

}  // namespace internal_r3x3_matrix

using internal_r3x3_matrix::KroneckerProduct;
using internal_r3x3_matrix::R3x3Matrix;

}  // namespace geometry
}  // namespace principia

#include "geometry/r3x3_matrix_body.hpp"
