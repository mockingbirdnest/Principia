
#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>

#include "geometry/r3_element.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace internal_r3x3_matrix {

using base::not_null;

// An |R3x3Matrix| is an element of the associative algebra of 3-by-3 matrices
// over ℝ, represented by |double|.
class R3x3Matrix {
 public:
  R3x3Matrix(R3Element<double> const& row_x,
             R3Element<double> const& row_y,
             R3Element<double> const& row_z);

  double Trace() const;
  R3x3Matrix Transpose() const;

  double operator()(int r, int c) const;

  R3x3Matrix& operator+=(R3x3Matrix const& right);
  R3x3Matrix& operator-=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(double right);
  R3x3Matrix& operator/=(double right);

  static R3x3Matrix Identity();

  void WriteToMessage(not_null<serialization::R3x3Matrix*> message) const;
  static R3x3Matrix ReadFromMessage(serialization::R3x3Matrix const& message);

 private:
  R3Element<double> row_x_;
  R3Element<double> row_y_;
  R3Element<double> row_z_;

  friend R3x3Matrix operator+(R3x3Matrix const& right);
  friend R3x3Matrix operator-(R3x3Matrix const& right);

  friend R3x3Matrix operator+(R3x3Matrix const& left,
                              R3x3Matrix const& right);
  friend R3x3Matrix operator-(R3x3Matrix const& left,
                              R3x3Matrix const& right);
  friend R3x3Matrix operator*(R3x3Matrix const& left,
                              R3x3Matrix const& right);

  template<typename Scalar>
  friend R3Element<Scalar> operator*(R3x3Matrix const& left,
                                     R3Element<Scalar> const& right);
  template<typename Scalar>
  friend R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                                     R3x3Matrix const& right);

  friend R3x3Matrix operator*(double left, R3x3Matrix const& right);
  friend R3x3Matrix operator*(R3x3Matrix const& left, double right);
  friend R3x3Matrix operator/(R3x3Matrix const& left, double right);

  friend bool operator==(R3x3Matrix const& left,
                         R3x3Matrix const& right);
  friend bool operator!=(R3x3Matrix const& left,
                         R3x3Matrix const& right);

  friend std::string DebugString(R3x3Matrix const& r3x3_matrix);
};

R3x3Matrix operator+(R3x3Matrix const& right);
R3x3Matrix operator-(R3x3Matrix const& right);

R3x3Matrix operator+(R3x3Matrix const& left,
                     R3x3Matrix const& right);
R3x3Matrix operator-(R3x3Matrix const& left,
                     R3x3Matrix const& right);
R3x3Matrix operator*(R3x3Matrix const& left,
                     R3x3Matrix const& right);

template<typename Scalar>
R3Element<Scalar> operator*(R3x3Matrix const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                            R3x3Matrix const& right);

R3x3Matrix operator*(double left, R3x3Matrix const& right);
R3x3Matrix operator*(R3x3Matrix const& left, double right);
R3x3Matrix operator/(R3x3Matrix const& left, double right);

bool operator==(R3x3Matrix const& left,
                R3x3Matrix const& right);
bool operator!=(R3x3Matrix const& left,
                R3x3Matrix const& right);

std::string DebugString(R3x3Matrix const& r3x3_matrix);

std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix const& r3x3_matrix);

}  // namespace internal_r3x3_matrix

using internal_r3x3_matrix::R3x3Matrix;

}  // namespace geometry
}  // namespace principia

#include "geometry/r3x3_matrix_body.hpp"
