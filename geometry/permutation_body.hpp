#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/Dimensionless.hpp"

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

// TODO(phl): Uncomment once orthogonal transformations are done.
/*template<typename FromFrame, typename ToFrame>
OrthogonalTransformation<FromFrame, ToFrame>
Permutation<FromFrame, ToFrame>::Forget() const {
  static const quantities::Dimentionless SqrtHalf(
      quantities::Dimentionless::Sqrt(quantities::Dimensionless(0.5)));
  static const R3Element<quantities::Dimensionless>[] = {
      R3Element<quantities::Dimensionless>(quantities::Dimensionless(0),
                                           quantities::Dimensionless(0),
                                           quantities::Dimensionless(0)),
      R3Element<Dimensionless>(
          Dimensionless(0.5), Dimensionless(0.5), Dimensionless(0.5)),
      R3Element<Dimensionless>(
          Dimensionless(0.5), Dimensionless(0.5), Dimensionless(0.5)),
      R3Element<Dimensionless>(Dimensionless(0), -SqrtHalf, SqrtHalf),
      R3Element<Dimensionless>(-SqrtHalf, SqrtHalf, Dimensionless(0)),
      R3Element<Dimensionless>(-SqrtHalf, Dimensionless(0), SqrtHalf)};
  static const Dimensionless[] quaternionRealParts = {
      Dimensionless(1),
      Dimensionless(-0.5),
      Dimensionless(0.5),
      Dimensionless(0),
      Dimensionless(0),
      Dimensionless(0)};
  const int i = 0x7 & (static_cast<int>(coordinate_permutation_) >> index);
  return OrthogonalTransformation<Scalar, FromFrame, ToFrame>(
      Determinant(),
      Rotation<A, B>(quaternionRealParts[i], quaternionImaginaryParts[i]));
}*/

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

}  // namespace geometry
}  // namespace principia
