#pragma once

#include "Grassmann.hpp"
#include "LinearMap.hpp"
#include "R3Element.hpp"
#include "Sign.hpp"
#include "..\Quantities\Dimensionless.hpp"

namespace Principia {
namespace Geometry {

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
Permutation<Scalar, FromFrame, ToFrame, Rank>::Permutation(
    const CoordinatePermutation coordinate_permutation)
    : coordinate_permutation_(coordinate_permutation) {}

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
Permutation<Scalar, FromFrame, ToFrame, Rank>::~Permutation() {}

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
MultiVector<Scalar, ToFrame, Rank> 
Permutation<Scalar, FromFrame, ToFrame, Rank>::ActOn(
    const MultiVector<Scalar, FromFrame, Rank>& right) {
}

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
Sign Permutation<Scalar, FromFrame, ToFrame, Rank>::Determinant() {
  return Sign(coordinate_permutation_ > 0); 
}

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
OrthogonalTransformation<Scalar, FromFrame, ToFrame, Rank> 
Permutation<Scalar, FromFrame, ToFrame, Rank>::Forget() {
  static const Quantities::Dimentionless SqrtHalf(Quantities::Dimentionless::Sqrt(Quantities::Dimensionless(0.5));
  static const R3Element<Quantities::Dimensionless>[] = {
      R3Element<Quantities::Dimensionless>(Quantities::Dimensionless(0), 
                                           Quantities::Dimensionless(0),
                                           Quantities::Dimensionless(0)),
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
  const int i = (0b111 & (static_cast<int>(coordinate_permutation_) >> index);
  return OrthogonalTransformation<Scalar, FromFrame, ToFrame, Rank>(
      Determinant(),
      Rotation<A, B>(quaternionRealParts[i], quaternionImaginaryParts[i]));
}

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>  
Permutation<Scalar, FromFrame, ToFrame, Rank>
Permutation<Scalar, FromFrame, ToFrame, Rank>::Identity() {
  return Permutation(XYZ);
}

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
R3Element<Scalar> operator*(
    const Permutation<Scalar, FromFrame, ToFrame, Rank>& left,
    const R3Element<Scalar>& right) {
  R3Element<Scalar> result;
  for (int coordinate = x; coordinate <= z; ++coordinate) {
    result[coordinate] = right[
        0b11 & 
        (static_cast<int>(left.coordinate_permutation_) >> (coordinate * 2))];
  }
  return result;
}

}  // namespace Geometry
}  // namespace Principia