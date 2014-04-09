#pragma once

namespace Principia {
namespace Geometry {
template<typename Scalar>
inline R3Element<Scalar>::R3Element(Scalar const& x,
                                    Scalar const& y,
                                    Scalar const& z) : x(x), y(y), z(z) {};


template<typename T>
inline R3Element<T> operator+(R3Element<T> const& right) {
  return R3Element<T>(+right.x, +right.y, +right.z);
}
template<typename T>
inline R3Element<T> operator-(R3Element<T> const& right) {
  return R3Element<T>(-right.x, -right.y, -right.z);
}

template<typename T>
inline R3Element<T> operator+(R3Element<T> const& left,
                              R3Element<T> const& right) {
  return R3Element<T>(left.x + right.x, left.y + right.y, left.z + right.z);
}
template<typename T>
inline R3Element<T> operator-(R3Element<T> const& left,
                              R3Element<T> const& right) {
  return R3Element<T>(left.x - right.x, left.y - right.y, left.z - right.z);
}

template<typename T>
inline R3Element<T> operator*(Quantities::Dimensionless const& left,
                              R3Element<T> const& right) {
  return R3Element<T>(left * right.x, left * right. y, left * right.z);
}
template<typename T>
inline R3Element<T> operator*(R3Element<T> const& left,
                              Quantities::Dimensionless const& right) {
  return R3Element<T>(left.x * right, left.y * right, left.z * right);
}
template<typename T>
inline R3Element<T> operator/(R3Element<T> const& left,
                              Quantities::Dimensionless const& right) {
  return R3Element<T>(left.x / right, left.y / right, left.z / right);
}

template<typename T, typename U>
inline R3Element<Quantities::Product<U, T>> operator*(
    U const& left,
    R3Element<T> const& right) {
  return R3Element<T>(left * right.x, left * right. y, left * right.z);
}
template<typename T, typename U>
inline R3Element<Quantities::Product<T, U>> operator*(R3Element<T> const& left,
                                                      U const& right) {
  return R3Element<T>(left.x * right, left.y * right, left.z * right);
}
template<typename T, typename U>
inline R3Element<Quantities::Quotient<T, U>> operator/(R3Element<T> const& left,
                                                       U const& right) {
  return R3Element<T>(left.x / right, left.y / right, left.z / right);
}

template<typename T>
inline void operator+=(R3Element<T>& left, R3Element<T> const& right) {
  left = left + right;
}
template<typename T>
inline void operator-=(R3Element<T>& left, R3Element<T> const& right) {
  right = left - right;
}

template<typename T>
inline void operator*=(R3Element<T>& left,
                       Quantities::Dimensionless const& right) {
  left = left * right;
}
template<typename T>
inline void operator/=(R3Element<T>& left,
                       Quantities::Dimensionless const& right) {
  left = left / right;
}

template<typename T, typename U>
inline R3Element<Quantities::Product<T, U>> Cross(R3Element<T> const& left,
                                                  R3Element<U> const& right) {
  return R3Element<Quantities::Product<T, U>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}
template<typename T, typename U>
inline Quantities::Product<T, U> Dot(R3Element<T> const& left,
                                     R3Element<U> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace Geometry
}  // namespace Principia
