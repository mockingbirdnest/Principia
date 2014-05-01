#pragma once

#include "Grassmann.hpp"
#include "LinearMap.hpp"
#include "R3Element.hpp"
#include "Sign.hpp"

namespace Principia {
namespace Geometry {

// A permutation of the coordinates. Obviously not coordinate-free, but
// practical. There are no precision losses when composing or applying
// permutations.
template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
class Permutation : public LinearMap<Scalar, FromFrame, ToFrame, Rank> {
  // The sign bit gives the sign of the permutation.
  static const int even = 0, odd = Int32.MinValue;
  // Three two-bit fields which indicate how each coordinate get mapped by the
  // permutation.
  static const int x = 0, y = 1, z = 2;
  // A three bit field used when using this enum to index arrays.
  static const int index = 6;
 public:
  enum CoordinatePermutation {
    XYZ = even + (x << x * 2) + (y << y * 2) + (z << z * 2) + (0 << index),
    YZX = even + (y << x * 2) + (z << y * 2) + (x << z * 2) + (1 << index),
    ZXY = even + (z << x * 2) + (x << y * 2) + (y << z * 2) + (2 << index),
    XZY = odd  + (x << x * 2) + (z << y * 2) + (y << z * 2) + (3 << index),
    ZYX = odd  + (z << x * 2) + (y << y * 2) + (x << z * 2) + (4 << index),
    YXZ = odd  + (y << x * 2) + (x << y * 2) + (z << z * 2) + (5 << index)
  }

  Permutation(const CoordinatePermutation coordinate_permutation);
  virtual ~Permutation();

  MultiVector<Scalar, ToFrame, Rank> ActOn(
      const MultiVector<Scalar, FromFrame, Rank>& right);

  Sign Determinant();

  OrthogonalTransformation<Scalar, FromFrame, ToFrame, Rank> Forget();

  static Permutation Identity();

 private:
  const CoordinatePermutation coordinate_permutation_;
};

template<typename Scalar, 
         typename FromFrame,
         typename ToFrame,
         unsigned int Rank>
R3Element<Scalar> operator*(
    const Permutation<Scalar, FromFrame, ToFrame, Rank>& left,
    const R3Element<Scalar>& right);

}  // namespace Geometry
}  // namespace Principia

#include "Permutation-body.hpp"