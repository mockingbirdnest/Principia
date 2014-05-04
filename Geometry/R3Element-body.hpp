#pragma once

namespace principia {
namespace geometry {
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
inline R3Element<Scalar> operator*(quantities::Dimensionless const& left,
                                   R3Element<Scalar> const& right) {
  return R3Element<Scalar>(left * right.x,
                           left * right.y,
                           left * right.z);
}
template<typename Scalar>
inline R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                                   quantities::Dimensionless const& right) {
  return R3Element<Scalar>(left.x * right,
                           left.y * right,
                           left.z * right);
}
template<typename Scalar>
inline R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                                   quantities::Dimensionless const& right) {
  return R3Element<Scalar>(left.x / right,
                           left.y / right,
                           left.z / right);
}

template<typename DLeft, typename Right>
inline R3Element<quantities::Product<quantities::Quantity<DLeft>,
                                     Right>> operator*(
    quantities::Quantity<DLeft> const& left,
    R3Element<Right> const& right) {
  return R3Element<quantities::Product<quantities::Quantity<DLeft>,
                                       Right>>(
      left * right.x,
      left * right. y,
      left * right.z);
}
template<typename LeftScalar, typename DRight>
inline R3Element<quantities::Product<LeftScalar,
                                     quantities::Quantity<DRight>>> operator*(
    R3Element<LeftScalar> const& left,
    quantities::Quantity<DRight> const& right) {
  return R3Element<quantities::Product<LeftScalar,
                                       quantities::Quantity<DRight>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}
template<typename LeftScalar, typename DRight>
inline R3Element<quantities::Quotient<LeftScalar,
                                      quantities::Quantity<DRight>>> operator/(
    R3Element<LeftScalar> const& left,
    quantities::Quantity<DRight> const& right) {
  return R3Element<quantities::Quotient<LeftScalar,
                                        quantities::Quantity<DRight>>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename Scalar>
inline void operator+=(R3Element<Scalar>& left,
                       R3Element<Scalar> const& right) {
  left = left + right;
}
template<typename Scalar>
inline void operator-=(R3Element<Scalar>& left,
                       R3Element<Scalar> const& right) {
  left = left - right;
}

template<typename Scalar>
inline void operator*=(R3Element<Scalar>& left,
                       quantities::Dimensionless const& right) {
  left = left * right;
}
template<typename Scalar>
inline void operator/=(R3Element<Scalar>& left,
                       quantities::Dimensionless const& right) {
  left = left / right;
}

template<typename LeftScalar, typename RightScalar>
inline R3Element<quantities::Product<LeftScalar, RightScalar>> Cross(
    R3Element<LeftScalar> const& left,
    R3Element<RightScalar> const& right) {
  return R3Element<quantities::Product<LeftScalar, RightScalar>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LeftScalar, typename RightScalar>
inline quantities::Product<LeftScalar, RightScalar> Dot(
    R3Element<LeftScalar> const& left,
    R3Element<RightScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace geometry
}  // namespace principia
