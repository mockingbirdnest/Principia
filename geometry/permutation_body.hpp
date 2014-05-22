#pragma once

#include <map>

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
Permutation<FromFrame, ToFrame>::Permutation(
    CoordinatePermutation const coordinate_permutation)
    : coordinate_permutation_(coordinate_permutation) {}

template<typename FromFrame, typename ToFrame>
inline Sign Permutation<FromFrame, ToFrame>::Determinant() const {
  return Sign(coordinate_permutation_);
}

template<typename FromFrame, typename ToFrame>
Permutation<ToFrame, FromFrame>
Permutation<FromFrame, ToFrame>::Inverse() const {
  static std::map<CoordinatePermutation, CoordinatePermutation> inverse = {
    {XYZ, XYZ},
    {YZX, ZXY},
    {ZXY, YZX},
    {XZY, XZY},
    {ZYX, ZYX},
    {YXZ, YXZ},
  };
  return Permutation(inverse[coordinate_permutation_]);
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>((*this)(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(
      Determinant() * (*this)(bivector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(Determinant() * trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> 
Permutation<FromFrame, ToFrame>::Forget() const {
  static const quantities::Dimensionless sqrt_half = quantities::Sqrt(0.5);
  static std::map<CoordinatePermutation, Quaternion> quaternion = {
    {XYZ, Quaternion(1, {0, 0, 0})},
    {YZX, Quaternion(0.5, {-0.5, -0.5, -0.5})},
    {ZXY, Quaternion(0.5, {0.5, 0.5, 0.5})},
    {XZY, Quaternion(0, {0, -sqrt_half, sqrt_half})},
    {ZYX, Quaternion(0, {-sqrt_half, 0, sqrt_half})},
    {YXZ, Quaternion(0, {-sqrt_half, sqrt_half, 0})},
  };
  return OrthogonalMap<FromFrame, ToFrame>(
      Determinant(),
      Rotation<FromFrame, ToFrame>(quaternion[coordinate_permutation_]));
}

template<typename FromFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> Permutation<FromFrame, ToFrame>::Identity() {
  return Permutation(XYZ);
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Permutation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  R3Element<Scalar> result;
  for (int coordinate = x; coordinate <= z; ++coordinate) {
    result[coordinate] = r3_element[
        0x3 &
        (static_cast<int>(coordinate_permutation_) >> (coordinate * 2))];
  }
  return result;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> operator*(
    Permutation<ThroughFrame, ToFrame> const& left,
    Permutation<FromFrame, ThroughFrame> const& right) {
  ///
}

}  // namespace geometry
}  // namespace principia
