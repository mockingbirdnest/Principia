#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/Quaternion.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Sign.hpp"
#include "Quantities/Dimensionless.hpp"
#include "Quantities/ElementaryFunctions.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Rotation<FromFrame, ToFrame>::Rotation(quantities::Angle const& angle,
                                       Vector<Scalar, FromFrame> const& axis) {
  quantities::Angle const half_angle = 0.5 * angle;
  quantities::Dimensionless const cos = Cos(half_angle);
  quantities::Dimensionless const sin = Sin(half_angle);
  R3Element<Scalar> const coordinates = axis.coordinates();
  Scalar const norm = coordinates.Norm();
  R3Element<quantities::Dimensionless> unit_axis = coordinates / norm;
  quaternion_ = Quaternion(cos, sin * unit_axis);
}

template<typename FromFrame, typename ToFrame>
Sign Rotation<FromFrame, ToFrame>::Determinant() const {
  return Sign(1);
}

template<typename FromFrame, typename ToFrame>
Rotation<ToFrame, FromFrame> Rotation<FromFrame, ToFrame>::Inverse() const {
  return Rotation(quaternion_.Conjugate());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Rotation<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>((*this)(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Rotation<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>((*this)(bivector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Rotation<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return trivector;
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::Identity() {
  return Rotation(Quaternion(1));
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame>::Rotation(Quaternion const& quaternion)
    : quaternion_(quaternion) {}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Rotation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  quantities::Dimensionless const& real_part = quaternion_.real_part();
  R3Element<quantities::Dimensionless> const& imaginary_part =
      quaternion_.imaginary_part();
  return r3_element + 2 * Cross(imaginary_part,
                                (Cross(imaginary_part, r3_element) +
                                    real_part * r3_element));
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right) {
  return Rotation<FromFrame, ToFrame>(left.quaternion_ * right.quaternion_);
}

}  // namespace geometry
}  // namespace principia
