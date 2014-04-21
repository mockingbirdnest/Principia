#pragma once

namespace Principia {
namespace Geometry {
template<typename Scalar>
inline R3Element<Scalar>::R3Element(Scalar const& x,
                                    Scalar const& y,
                                    Scalar const& z) : x(x), y(y), z(z) {};


template<typename Scalar>
inline R3Element<Scalar> operator+(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(+right.x, +right.y, +right.z);
}
template<typename Scalar>
inline R3Element<Scalar> operator-(R3Element<Scalar> const& right) {
  return R3Element<Scalar>(-right.x, -right.y, -right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator+(
    R3Element<Scalar> const& left,
    R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x + right.x,
                           left.y + right.y,
                           left.z + right.z);
}
template<typename Scalar>
inline R3Element<Scalar> operator-(
    R3Element<Scalar> const& left,
    R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left.x - right.x,
                           left.y - right.y,
                           left.z - right.z);
}

template<typename Scalar>
inline R3Element<Scalar> operator*(Quantities::Dimensionless const& left,
                                   R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left * right.x,
                           left * right.y,
                           left * right.z);
}
template<typename Scalar>
inline R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                                   Quantities::Dimensionless const& right) {
  return R3Element<Scalar>(left.x * right,
                           left.y * right,
                           left.z * right);
}
template<typename Scalar>
inline R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                                   Quantities::Dimensionless const& right) {
  return R3Element<Scalar>(left.x / right,
                           left.y / right,
                           left.z / right);
}

template<typename DLeft, typename DRight>
inline R3Element<Quantities::Product<Quantities::Quantity<DLeft>,
                                     Quantities::Quantity<DRight>>> operator*(
    Quantities::Quantity<DRight> const& left,
    R3Element<Quantities::Quantity<DRight>> const& right) {
  return R3Element<Quantities::Product<Quantities::Quantity<DLeft>,
                                       Quantities::Quantity<DRight>>>(
      left * right.x,
      left * right. y,
      left * right.z);
}
template<typename DLeft, typename DRight>
inline R3Element<Quantities::Product<Quantities::Quantity<DLeft>,
                                     Quantities::Quantity<DRight>>> operator*(
    R3Element<Quantities::Quantity<DLeft>> const& left,
    Quantities::Quantity<DRight> const& right) {
  return R3Element<Quantities::Product<Quantities::Quantity<DLeft>,
                                       Quantities::Quantity<DRight>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}
template<typename DLeft, typename DRight>
inline R3Element<Quantities::Quotient<Quantities::Quantity<DLeft>,
                                      Quantities::Quantity<DRight>>> operator/(
    R3Element<Quantities::Quantity<DLeft>> const& left,
    Quantities::Quantity<DRight> const& right) {
  return R3Element<Quantities::Quotient<Quantities::Quantity<DLeft>,
                                        Quantities::Quantity<DRight>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename D>
inline void operator+=(R3Element<Quantities::Quantity<D>>& left,
                       R3Element<Quantities::Quantity<D>> const& right) {
  left = left + right;
}
template<typename D>
inline void operator-=(R3Element<Quantities::Quantity<D>>& left,
                       R3Element<Quantities::Quantity<D>> const& right) {
  left = left - right;
}

template<typename D>
inline void operator*=(R3Element<Quantities::Quantity<D>>& left,
                       Quantities::Dimensionless const& right) {
  left = left * right;
}
template<typename D>
inline void operator/=(R3Element<Quantities::Quantity<D>>& left,
                       Quantities::Dimensionless const& right) {
  left = left / right;
}

template<typename Left, typename Right>
inline R3Element<Quantities::Product<Left, Right>> Cross(
    R3Element<Left> const& left,
    R3Element<Right> const& right) {
  return R3Element<Quantities::Product<Left, Right>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename DLeft, typename DRight>
inline Quantities::Product<Quantities::Quantity<DLeft>,
                           Quantities::Quantity<DRight>> Dot(
    R3Element<Quantities::Quantity<DLeft>> const& left,
    R3Element<Quantities::Quantity<DRight>> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace Geometry
}  // namespace Principia
