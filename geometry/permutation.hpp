#pragma once

#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "gtest/gtest_prod.h"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class OrthogonalMap;

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
  ~Permutation() override = default;

  Sign Determinant() const override;

  Permutation<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename T>
  typename base::Mappable<Permutation, T>::type operator()(T const& t) const; 

  OrthogonalMap<FromFrame, ToFrame> Forget() const;

  static Permutation Identity();

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  CoordinatePermutation coordinate_permutation_;

  template<typename From, typename Through, typename To>
  friend Permutation<From, To> operator*(
      Permutation<Through, To> const& left,
      Permutation<From, Through> const& right);

  // As much as I dislike FRIEND_TEST(), it seems like the most convenient way
  // to access the above operator.
  FRIEND_TEST(PermutationTest, Identity);
  FRIEND_TEST(PermutationTest, XYZ);
  FRIEND_TEST(PermutationTest, YZX);
  FRIEND_TEST(PermutationTest, ZXY);
  FRIEND_TEST(PermutationTest, XZY);
  FRIEND_TEST(PermutationTest, ZYX);
  FRIEND_TEST(PermutationTest, YXZ);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> operator*(
    Permutation<ThroughFrame, ToFrame> const& left,
    Permutation<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/permutation_body.hpp"
