
#pragma once

#include "geometry/r3x3_matrix.hpp"

#include <string>
#include <utility>

#include "base/macros.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_r3x3_matrix {

using quantities::SIUnit;

template<typename Scalar>
R3x3Matrix<Scalar>::R3x3Matrix(R3Element<Scalar> const& row_x,
                               R3Element<Scalar> const& row_y,
                               R3Element<Scalar> const& row_z)
    : rows_({row_x, row_y, row_z}) {}

template<typename Scalar>
R3x3Matrix<Scalar> R3x3Matrix<Scalar>::Diagonal(
    R3Element<Scalar> const& diagonal) {
  return {
      {diagonal.x, {}, {}},
      {{}, diagonal.y, {}},
      {{}, {}, diagonal.z},
  };
}

template<typename Scalar>
Scalar R3x3Matrix<Scalar>::Trace() const {
  return rows_[X].x + rows_[Y].y + rows_[Z].z;
}

template<typename Scalar>
Cube<Scalar> R3x3Matrix<Scalar>::Determinant() const {
  return rows_[X].x * rows_[Y].y * rows_[Z].z +
         rows_[X].y * rows_[Y].z * rows_[Z].x +
         rows_[X].z * rows_[Y].x * rows_[Z].y -
         rows_[X].z * rows_[Y].y * rows_[Z].x -
         rows_[X].y * rows_[Y].x * rows_[Z].z -
         rows_[X].x * rows_[Y].z * rows_[Z].y;
}

template<typename Scalar>
R3x3Matrix<Scalar> R3x3Matrix<Scalar>::Transpose() const {
  return R3x3Matrix({rows_[X].x, rows_[Y].x, rows_[Z].x},
                    {rows_[X].y, rows_[Y].y, rows_[Z].y},
                    {rows_[X].z, rows_[Y].z, rows_[Z].z});
}

template<typename Scalar>
R3Element<Scalar> const& R3x3Matrix<Scalar>::row_x() const {
  return rows_[X];
}

template<typename Scalar>
R3Element<Scalar> const& R3x3Matrix<Scalar>::row_y() const {
  return rows_[Y];
}

template<typename Scalar>
R3Element<Scalar> const& R3x3Matrix<Scalar>::row_z() const {
  return rows_[Z];
}

template<typename Scalar>
FORCE_INLINE(inline) Scalar R3x3Matrix<Scalar>::operator()(
    int const r, int const c) const {
  switch (r) {
    case 0:
    case 1:
    case 2:
      return rows_[r][c];
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE
                  << " indices = {" << r << ", " << c << "}";
      base::noreturn();
  }
}

template<typename Scalar>
Scalar& R3x3Matrix<Scalar>::operator()(int r, int c) {
  switch (r) {
    case 0:
    case 1:
    case 2:
      return rows_[r][c];
    default:
      DLOG(FATAL) << FUNCTION_SIGNATURE
                  << " indices = {" << r << ", " << c << "}";
      base::noreturn();
  }
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator+=(
    R3x3Matrix const& right) {
  return *this = *this + right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator-=(
    R3x3Matrix const& right) {
  return *this = *this - right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator*=(
    R3x3Matrix const& right) {
  return *this = *this * right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator*=(double const right) {
  return *this = *this * right;
}

template<typename Scalar>
R3x3Matrix<Scalar>& R3x3Matrix<Scalar>::operator/=(double const right) {
  return *this = *this / right;
}

template<typename Scalar>
template<typename RScalar>
R3Element<Quotient<RScalar, Scalar>> R3x3Matrix<Scalar>::Solve(
    R3Element<RScalar> const& rhs) const {
  R3x3Matrix<double> l;
  R3x3Matrix<Scalar> u;

  // Doolittle's method.
  for (int k = 0; k < 3; ++k) {
    // Partial pivoting.
    int r = -1;
    Scalar max;
    for (int i = k; i < 3; ++i) {
      if (Abs((*this)(i, k)) > max) {
        r = i;
      }
    }
    //swap rows k and r.
    //swap(rhs[k], rhs[r]);

    for (int j = k; j < 3; ++j) {
      Scalar u_kj = (*this)(k, j);
      for (int i = 0; i < k - 1; ++i) {
        u_kj -= l(k, i) * u(i, j)
      }
      u(k, j) = u_kj;
    }
    for (int i = k + 1; i < 3; ++i) {
      Scalar l_ik = (*this)(i, k);
      for (int j = 0; j < k - 1; ++j) {
        l_ik -= l(i, j) * u (j, k)
      }
      l(i, k) = l_ik / u(k, k)
    }
  }
}

template<typename Scalar>
template<typename S, typename>
R3x3Matrix<S> R3x3Matrix<Scalar>::Identity() {
  return R3x3Matrix<S>({1, 0, 0},
                       {0, 1, 0},
                       {0, 0, 1});
}

template<typename Scalar>
void R3x3Matrix<Scalar>::WriteToMessage(
    not_null<serialization::R3x3Matrix*> const message) const {
  rows_[X].WriteToMessage(message->mutable_row_x());
  rows_[Y].WriteToMessage(message->mutable_row_y());
  rows_[Z].WriteToMessage(message->mutable_row_z());
}

template<typename Scalar>
R3x3Matrix<Scalar> R3x3Matrix<Scalar>::ReadFromMessage(
    serialization::R3x3Matrix const& message) {
  return R3x3Matrix(R3Element<Scalar>::ReadFromMessage(message.row_x()),
                    R3Element<Scalar>::ReadFromMessage(message.row_y()),
                    R3Element<Scalar>::ReadFromMessage(message.row_z()));
}

template<typename Scalar>
R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& right) {
  constexpr auto X = R3x3Matrix<Scalar>::X;
  constexpr auto Y = R3x3Matrix<Scalar>::Y;
  constexpr auto Z = R3x3Matrix<Scalar>::Z;
  return R3x3Matrix<Scalar>(+right.rows_[X], +right.rows_[Y], +right.rows_[Z]);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& right) {
  constexpr auto X = R3x3Matrix<Scalar>::X;
  constexpr auto Y = R3x3Matrix<Scalar>::Y;
  constexpr auto Z = R3x3Matrix<Scalar>::Z;
  return R3x3Matrix<Scalar>(-right.rows_[X], -right.rows_[Y], -right.rows_[Z]);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator+(R3x3Matrix<Scalar> const& left,
                             R3x3Matrix<Scalar> const& right) {
  constexpr auto X = R3x3Matrix<Scalar>::X;
  constexpr auto Y = R3x3Matrix<Scalar>::Y;
  constexpr auto Z = R3x3Matrix<Scalar>::Z;
  return R3x3Matrix<Scalar>(left.rows_[X] + right.rows_[X],
                            left.rows_[Y] + right.rows_[Y],
                            left.rows_[Z] + right.rows_[Z]);
}

template<typename Scalar>
R3x3Matrix<Scalar> operator-(R3x3Matrix<Scalar> const& left,
                             R3x3Matrix<Scalar> const& right) {
  constexpr auto X = R3x3Matrix<Scalar>::X;
  constexpr auto Y = R3x3Matrix<Scalar>::Y;
  constexpr auto Z = R3x3Matrix<Scalar>::Z;
  return R3x3Matrix<Scalar>(left.rows_[X] - right.rows_[X],
                            left.rows_[Y] - right.rows_[Y],
                            left.rows_[Z] - right.rows_[Z]);
}

template<typename LScalar, typename RScalar>
R3x3Matrix<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    R3x3Matrix<RScalar> const& right) {
  constexpr auto X = R3x3Matrix<LScalar>::X;
  constexpr auto Y = R3x3Matrix<LScalar>::Y;
  constexpr auto Z = R3x3Matrix<LScalar>::Z;
  R3x3Matrix<RScalar> const t_right = right.Transpose();
  return R3x3Matrix<Product<LScalar, RScalar>>(
             {Dot(left.rows_[X], t_right.rows_[X]),
              Dot(left.rows_[X], t_right.rows_[Y]),
              Dot(left.rows_[X], t_right.rows_[Z])},
             {Dot(left.rows_[Y], t_right.rows_[X]),
              Dot(left.rows_[Y], t_right.rows_[Y]),
              Dot(left.rows_[Y], t_right.rows_[Z])},
             {Dot(left.rows_[Z], t_right.rows_[X]),
              Dot(left.rows_[Z], t_right.rows_[Y]),
              Dot(left.rows_[Z], t_right.rows_[Z])});
}

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    R3Element<RScalar> const& right) {
  constexpr auto X = R3x3Matrix<LScalar>::X;
  constexpr auto Y = R3x3Matrix<LScalar>::Y;
  constexpr auto Z = R3x3Matrix<LScalar>::Z;
  return R3Element<Product<LScalar, RScalar>>({Dot(left.rows_[X], right),
                                               Dot(left.rows_[Y], right),
                                               Dot(left.rows_[Z], right)});
}

template<typename LScalar, typename RScalar>
R3Element<Product<LScalar, RScalar>> operator*(
    R3Element<LScalar> const& left,
    R3x3Matrix<RScalar> const& right) {
  R3x3Matrix<RScalar> const t_right = right.Transpose();
  constexpr auto X = R3x3Matrix<LScalar>::X;
  constexpr auto Y = R3x3Matrix<LScalar>::Y;
  constexpr auto Z = R3x3Matrix<LScalar>::Z;
  return R3Element<Product<LScalar, RScalar>>({Dot(left, t_right.rows_[X]),
                                               Dot(left, t_right.rows_[Y]),
                                               Dot(left, t_right.rows_[Z])});
}


template<typename LScalar, typename RScalar, typename>
R3x3Matrix<Product<LScalar, RScalar>> operator*(
    LScalar const& left,
    R3x3Matrix<RScalar> const& right) {
  constexpr auto X = R3x3Matrix<RScalar>::X;
  constexpr auto Y = R3x3Matrix<RScalar>::Y;
  constexpr auto Z = R3x3Matrix<RScalar>::Z;
  return R3x3Matrix<Product<LScalar, RScalar>>(left * right.rows_[X],
                                               left * right.rows_[Y],
                                               left * right.rows_[Z]);
}

template<typename LScalar, typename RScalar, typename>
R3x3Matrix<Product<LScalar, RScalar>> operator*(
    R3x3Matrix<LScalar> const& left,
    RScalar const& right) {
  constexpr auto X = R3x3Matrix<LScalar>::X;
  constexpr auto Y = R3x3Matrix<LScalar>::Y;
  constexpr auto Z = R3x3Matrix<LScalar>::Z;
  return R3x3Matrix<Product<LScalar, RScalar>>(left.rows_[X] * right,
                                               left.rows_[Y] * right,
                                               left.rows_[Z] * right);
}

template<typename LScalar, typename RScalar, typename>
R3x3Matrix<Quotient<LScalar, RScalar>> operator/(
    R3x3Matrix<LScalar> const& left,
    RScalar const& right) {
  constexpr auto X = R3x3Matrix<LScalar>::X;
  constexpr auto Y = R3x3Matrix<LScalar>::Y;
  constexpr auto Z = R3x3Matrix<LScalar>::Z;
  return R3x3Matrix<Quotient<LScalar, RScalar>>(left.rows_[X] / right,
                                                left.rows_[Y] / right,
                                                left.rows_[Z] / right);
}

template<typename LScalar, typename RScalar>
R3x3Matrix<Product<LScalar, RScalar>> KroneckerProduct(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3x3Matrix<Product<LScalar, RScalar>>(left.x * right,
                                               left.y * right,
                                               left.z * right);
}

template<typename Scalar>
bool operator==(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right) {
  return left.rows_ == right.rows_;
}

template<typename Scalar>
bool operator!=(R3x3Matrix<Scalar> const& left,
                R3x3Matrix<Scalar> const& right) {
  return left.rows_ != right.rows_;
}

template<typename Scalar>
std::string DebugString(R3x3Matrix<Scalar> const& r3x3_matrix) {
  constexpr auto X = R3x3Matrix<Scalar>::X;
  constexpr auto Y = R3x3Matrix<Scalar>::Y;
  constexpr auto Z = R3x3Matrix<Scalar>::Z;
  std::string result = "{";
  result += DebugString(r3x3_matrix.rows_[X]);
  result += ", ";
  result += DebugString(r3x3_matrix.rows_[Y]);
  result += ", ";
  result += DebugString(r3x3_matrix.rows_[Z]);
  result +="}";
  return result;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3x3Matrix<Scalar> const& r3x3_matrix) {
  out << DebugString(r3x3_matrix);
  return out;
}

}  // namespace internal_r3x3_matrix
}  // namespace geometry
}  // namespace principia
