#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Sign.hpp"
#include "quantities/Dimensionless.hpp"

namespace principia {
namespace geometry {

template<typename Scalar, typename FromFrame, typename ToFrame>
Permutation<Scalar, FromFrame, ToFrame>::Permutation(
    CoordinatePermutation const coordinate_permutation)
    : coordinate_permutation_(coordinate_permutation) {}

template<typename Scalar, typename FromFrame, typename ToFrame>
Vector<Scalar, ToFrame> Permutation<Scalar, FromFrame, ToFrame>::ActOn(
    Vector<Scalar, FromFrame> const& vector) const {
  return *this * vector.coordinates;
}

template<typename Scalar, typename FromFrame, typename ToFrame>
Bivector<Scalar, ToFrame> Permutation<Scalar, FromFrame, ToFrame>::ActOn(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Determinant() * (*this * bivector.coordinates);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
Trivector<Scalar, ToFrame> Permutation<Scalar, FromFrame, ToFrame>::ActOn(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Determinant() * trivector.coordinates;
}

template<typename Scalar, typename FromFrame, typename ToFrame>
Sign Permutation<Scalar, FromFrame, ToFrame>::Determinant() const {
  return Sign(coordinate_permutation_);
}

/*template<typename Scalar, typename FromFrame, typename ToFrame>
OrthogonalTransformation<Scalar, FromFrame, ToFrame> 
Permutation<Scalar, FromFrame, ToFrame>::Forget() const {
  static const quantities::Dimentionless SqrtHalf(quantities::Dimentionless::Sqrt(quantities::Dimensionless(0.5)));
  static const R3Element<quantities::Dimensionless>[] = {
      R3Element<quantities::Dimensionless>(quantities::Dimensionless(0), 
                                           quantities::Dimensionless(0),
                                           quantities::Dimensionless(0)),
      R3Element<Dimensionless>(Dimensionless(0.5), Dimensionless(0.5), Dimensionless(0.5)),
      R3Element<Dimensionless>(Dimensionless(0.5), Dimensionless(0.5), Dimensionless(0.5)),
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

template<typename Scalar, typename FromFrame, typename ToFrame>
Permutation<Scalar, FromFrame, ToFrame>
Permutation<Scalar, FromFrame, ToFrame>::Identity() {
  return Permutation(XYZ);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
R3Element<Scalar> operator*(
    Permutation<Scalar, FromFrame, ToFrame> const& left,
    R3Element<Scalar> const& right) {
  typedef Permutation<Scalar, FromFrame, ToFrame> P;
  R3Element<Scalar> result;
  for (int coordinate = P::x; coordinate <= P::z; ++coordinate) {
    result[coordinate] = right[
        0x3 & 
        (static_cast<int>(left.coordinate_permutation_) >> (coordinate * 2))];
  }
  return result;
}

}  // namespace geometry
}  // namespace principia
