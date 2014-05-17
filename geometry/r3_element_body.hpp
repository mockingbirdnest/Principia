#pragma once

#include <assert.h>

namespace principia {
namespace geometry {

template<typename Scalar>
inline R3Element<Scalar>::R3Element(Scalar const& x,
                                    Scalar const& y,
                                    Scalar const& z) : x(x), y(y), z(z) {}

template<typename Scalar>
inline Scalar& R3Element<Scalar>::operator[](int const index) {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      // TODO(phl): Use glog?
      assert(false);
      return x;
  }
}

template<typename Scalar>
inline Scalar const& R3Element<Scalar>::operator[](int const index) const {
  switch (index) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      // TODO(phl): Use glog?
      assert(false);
      return x;
  }
}

template<typename Scalar>
inline Scalar R3Element<Scalar>::Norm() const {
  return Sqrt(Dot(*this, *this));
}

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

template<typename LDimension, typename RScalar>
inline R3Element<quantities::Product<quantities::Quantity<LDimension>, RScalar>>
operator*(quantities::Quantity<LDimension> const& left,
          R3Element<RScalar> const& right) {
  return R3Element<quantities::Product<quantities::Quantity<LDimension>,
                                       RScalar>>(
      left * right.x,
      left * right.y,
      left * right.z);
}

template<typename LScalar, typename RDimension>
inline R3Element<quantities::Product<LScalar, quantities::Quantity<RDimension>>>
operator*(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3Element<quantities::Product<LScalar,
                                       quantities::Quantity<RDimension>>>(
      left.x * right,
      left.y * right,
      left.z * right);
}

template<typename LScalar, typename RDimension>
inline R3Element<quantities::Quotient<LScalar,
                                      quantities::Quantity<RDimension>>>
operator/(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right) {
  return R3Element<quantities::Quotient<LScalar,
                                        quantities::Quantity<RDimension>>>(
      left.x / right,
      left.y / right,
      left.z / right);
}

template<typename Scalar>
bool operator==(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template<typename Scalar>
bool operator!=(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right) {
  return left.x != right.x || left.y != right.y || left.z != right.z;
}

template<typename Scalar>
inline void operator+=(R3Element<Scalar>& left,  // NOLINT(runtime/references)
                       R3Element<Scalar> const& right) {
  left = left + right;
}

template<typename Scalar>
inline void operator-=(R3Element<Scalar>& left,  // NOLINT(runtime/references)
                       R3Element<Scalar> const& right) {
  left = left - right;
}

template<typename Scalar>
inline void operator*=(R3Element<Scalar>& left,  // NOLINT(runtime/references)
                       quantities::Dimensionless const& right) {
  left = left * right;
}

template<typename Scalar>
inline void operator/=(R3Element<Scalar>& left,  // NOLINT(runtime/references)
                       quantities::Dimensionless const& right) {
  left = left / right;
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3Element<Scalar> const& r3_element) {
  out << "{" << r3_element.x << ", "
             << r3_element.y << ", "
             << r3_element.z << "}";
  return out;
}

template<typename LScalar, typename RScalar>
inline R3Element<quantities::Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return R3Element<quantities::Product<LScalar, RScalar>>(
      left.y * right.z - left.z * right.y,
      left.z * right.x - left.x * right.z,
      left.x * right.y - left.y * right.x);
}

template<typename LScalar, typename RScalar>
inline quantities::Product<LScalar, RScalar> Dot(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

}  // namespace geometry
}  // namespace principia
