#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {

// A permutation of the coordinates. Obviously not coordinate-free, but
// practical.  There are no precision losses when composing or applying
// permutations.
template<typename FromFrame, typename ToFrame>
class Permutation : public LinearMap<FromFrame, ToFrame> {
  // These constants are used in the definition of type CoordinatePermutation.
  // The sign bit gives the sign of the permutation.
  static int const even = 0, odd = 0x80000000;
  // Three two-bit fields which indicate how each coordinate get mapped by the
  // permutation.
  static int const x = 0, y = 1, z = 2;
  // A three bit field used when using this enum to index arrays.
  static int const index = 6;

 public:
  enum CoordinatePermutation {
    XYZ = even + (x << x * 2) + (y << y * 2) + (z << z * 2) + (0 << index),
    YZX = even + (y << x * 2) + (z << y * 2) + (x << z * 2) + (1 << index),
    ZXY = even + (z << x * 2) + (x << y * 2) + (y << z * 2) + (2 << index),
    XZY = odd  + (x << x * 2) + (z << y * 2) + (y << z * 2) + (3 << index),
    ZYX = odd  + (z << x * 2) + (y << y * 2) + (x << z * 2) + (4 << index),
    YXZ = odd  + (y << x * 2) + (x << y * 2) + (z << z * 2) + (5 << index)
  };

  explicit Permutation(CoordinatePermutation const coordinate_permutation);
  virtual ~Permutation() = default;

  Sign Determinant() const override;

  // TODO(phl): Inverse, composition.

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  // TODO(phl): Uncomment once orthogonal transformations are done.
  // OrthogonalTransformation<Scalar, FromFrame, ToFrame> Forget() const;

  static Permutation Identity();

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  CoordinatePermutation const coordinate_permutation_;
  friend class PermutationTests;
};

}  // namespace geometry
}  // namespace principia

#include "geometry/permutation_body.hpp"
