#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Sign.hpp"
#include "Quantities/Dimensionless.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame>::Rotation(
    quantities::Dimensionless const& real_part,
    R3Element<quantities::Dimensionless> const& imaginary_part)
    : real_part_(real_part),
      imaginary_part_(imaginary_part) {}

template<typename FromFrame, typename ToFrame>
Sign Rotation<FromFrame, ToFrame>::Determinant() const {
  return Sign(1);
}

template<typename FromFrame, typename ToFrame>
Rotation<ToFrame, FromFrame> Rotation<FromFrame, ToFrame>::Inverse() const {
  return Rotation(real_part_, -imaginary_part_);
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
  return Rotation(1, R3Element<quantities::Dimensionless>(0, 0, 0));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Rotation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  // TODO(egg): Optimise.
  // TODO(phl): WTF
  return r3_element;

  //return Maps.Compose<B, A, B>(
  //  Maps.Compose<A, A, B>(left, new Rotation<A, A>((Scalar)0, right)),
  //  left.Inverse()).imaginary_part_;
}

template<typename FromFrame, typename ToFrame>
void Rotation<FromFrame, ToFrame>::QuaternionMultiplication(
    quantities::Dimensionless const& left_real_part,
    R3Element<quantities::Dimensionless> const& left_imaginary_part,
    quantities::Dimensionless const& right_real_part,
    R3Element<quantities::Dimensionless> const& right_imaginary_part,
    quantities::Dimensionless* result_real_part,
    R3Element<quantities::Dimensionless>* result_imaginary_part) {
  *result_real_part = left_real_part_ * right_real_part_ -
      left_imaginary_part_.Dot(right_imaginary_part_);
  *result_imaginary_part =
      left.real_part_ * right.imaginary_part_ +
          right.real_part_ * left.imaginary_part_ +
          left.imaginary_part_.Cross(right.imaginary_part_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right) {
  return Rotation(
      left.real_part_ * right.real_part_ - 
          left.imaginary_part_.Dot(right.imaginary_part_),
      left.real_part_ * right.imaginary_part_ + 
          right.real_part_ * left.imaginary_part_ + 
          left.imaginary_part_.Cross(right.imaginary_part_));
}

}  // namespace geometry
}  // namespace principia
