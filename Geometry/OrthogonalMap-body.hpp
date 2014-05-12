#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/OrthogonalMap.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Sign.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
OrthogonalMap<FromFrame, ToFrame>::OrthogonalMap(
    const Sign& determinant,
    Rotation<FromFrame, ToFrame> const& rotation)
    : determinant_(determinant),
      rotation_(rotation) {}

template<typename FromFrame, typename ToFrame>
Sign OrthogonalMap<FromFrame, ToFrame>::Determinant() const {
  return determinant_;
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<ToFrame, FromFrame> OrthogonalMap<FromFrame, ToFrame>::Inverse() const {
  return OrthogonalMap(determinant_, rotation_.Inverse());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>((*this)(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(this->rotation_(bivector));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return determinant_ * trivector;
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> OrthogonalMap<FromFrame, ToFrame>::Identity() {
  return OrthogonalMap(Sign(1), Rotation<FromFrame, ToFrame>::Identity());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> OrthogonalMap<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  return rotation_(determinant_ * r3_element);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right) {
  return OrthogonalMap(left.determinant_ * right.determinant_,
                       left.rotation_ * right.rotation_);
}

}  // namespace geometry
}  // namespace principia
