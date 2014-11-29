#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {

namespace {

Quaternion ToQuaternion(R3x3Matrix const& matrix) {
  double const t = matrix.Trace();
  CHECK_LE(0, t);
  double const r = sqrt(1.0 + t);
  double const s = 0.5 / r;
  double const real_part = 0.5 * r;
  R3Element<double> imaginary_part;
  imaginary_part.x = (matrix[{2, 1}] - matrix[{1, 2}]) * s;
  imaginary_part.y = (matrix[{0, 2}] - matrix[{2, 0}]) * s;
  imaginary_part.z = (matrix[{1, 0}] - matrix[{0, 1}]) * s;
  return Quaternion(real_part, imaginary_part);
}

}  // namespace

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame>::Rotation() : quaternion_(Quaternion(1)) {}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame>::Rotation(Quaternion const& quaternion)
    : quaternion_(quaternion) {}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame>::Rotation(R3x3Matrix const& matrix)
  : Rotation(ToQuaternion(matrix)) {}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Rotation<FromFrame, ToFrame>::Rotation(
    quantities::Angle const& angle,
    Bivector<Scalar, FromFrame> const& axis) {
  quantities::Angle const half_angle = 0.5 * angle;
  double const cos = Cos(half_angle);
  double const sin = Sin(half_angle);
  R3Element<Scalar> const coordinates = axis.coordinates();
  Scalar const norm = coordinates.Norm();
  R3Element<double> const unit_axis = coordinates / norm;
  quaternion_ = Quaternion(cos, sin * unit_axis);
}

template<typename FromFrame, typename ToFrame>
Sign Rotation<FromFrame, ToFrame>::Determinant() const {
  return Sign(1);
}

template<typename FromFrame, typename ToFrame>
Rotation<ToFrame, FromFrame> Rotation<FromFrame, ToFrame>::Inverse() const {
  // Because |quaternion_| has norm 1, its inverse is just its conjugate.
  return Rotation<ToFrame, FromFrame>(quaternion_.Conjugate());
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
OrthogonalMap<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::Forget() const {
  return OrthogonalMap<FromFrame, ToFrame>(Sign(1), *this);
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> Rotation<FromFrame, ToFrame>::Identity() {
  return Rotation(Quaternion(1));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Rotation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  double const real_part = quaternion_.real_part();
  R3Element<double> const& imaginary_part = quaternion_.imaginary_part();
  return r3_element + 2 * Cross(imaginary_part,
                                Cross(imaginary_part, r3_element) +
                                    real_part * r3_element);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right) {
  return Rotation<FromFrame, ToFrame>(left.quaternion_ * right.quaternion_);
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Rotation<FromFrame, ToFrame> const& rotation) {
  return out << rotation.quaternion_;
}

}  // namespace geometry
}  // namespace principia
