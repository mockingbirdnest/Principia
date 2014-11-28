#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>

#include "geometry/r3_element.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {

// An |R3x3Matrix| is an element of the associative algebra of 3-by-3 matrices
// over ℝ, represented by |double|.
class R3x3Matrix {
 public:
  // The identity matrix.
  R3x3Matrix();
  R3x3Matrix(R3Element<double> const& column_x,
             R3Element<double> const& column_y,
             R3Element<double> const& column_z);

  R3x3Matrix& operator+=(R3x3Matrix const& right);
  R3x3Matrix& operator-=(R3x3Matrix const& right);
  R3x3Matrix& operator*=(double const right);
  R3x3Matrix& operator/=(double const right);

 private:
  R3Element<double> column_x_;
  R3Element<double> column_y_;
  R3Element<double> column_z_;

  friend R3x3Matrix operator+(R3x3Matrix const& right);
  friend R3x3Matrix operator-(R3x3Matrix const& right);

  friend R3x3Matrix operator+(R3x3Matrix const& left,
                              R3x3Matrix const& right);
  friend R3x3Matrix operator-(R3x3Matrix const& left,
                              R3x3Matrix const& right);

  friend R3x3Matrix operator*(double const left,
                              R3x3Matrix const& right);
  friend R3x3Matrix operator*(R3x3Matrix const& left,
                              double const right);
  friend R3x3Matrix operator/(R3x3Matrix const& left,
                              double const right);

  friend bool operator==(R3x3Matrix const& left,
                         R3x3Matrix const& right);
  friend bool operator!=(R3x3Matrix const& left,
                         R3x3Matrix const& right);

  friend std::string DebugString(R3x3Matrix const& r3_element);
};

R3x3Matrix operator+(R3x3Matrix const& right);
R3x3Matrix operator-(R3x3Matrix const& right);

R3x3Matrix operator+(R3x3Matrix const& left,
                     R3x3Matrix const& right);
R3x3Matrix operator-(R3x3Matrix const& left,
                     R3x3Matrix const& right);

R3x3Matrix operator*(double const left,
                     R3x3Matrix const& right);
R3x3Matrix operator*(R3x3Matrix const& left,
                     double const right);
R3x3Matrix operator/(R3x3Matrix const& left,
                     double const right);

bool operator==(R3x3Matrix const& left,
                R3x3Matrix const& right);
bool operator!=(R3x3Matrix const& left,
                R3x3Matrix const& right);

std::string DebugString(R3x3Matrix const& r3_element);

std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix const& r3_element);

}  // namespace geometry
}  // namespace principia

#include "geometry/r3x3_matrix_body.hpp"
