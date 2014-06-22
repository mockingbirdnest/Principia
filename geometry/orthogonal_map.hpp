#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {

// An orthogonal map between the inner product spaces |FromFrame| and
// |ToFrame|, as well as the induced maps on the exterior algebra.
// The orthogonal map is modeled as a rotoinversion.
template<typename FromFrame, typename ToFrame>
class OrthogonalMap : public LinearMap<FromFrame, ToFrame> {
 public:
  OrthogonalMap();
  virtual ~OrthogonalMap() = default;

  Sign Determinant() const override;

  OrthogonalMap<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  static OrthogonalMap Identity();

 private:
  OrthogonalMap(Sign const& determinant,
                Rotation<FromFrame, ToFrame> const& rotation);

  Sign determinant_;
  Rotation<FromFrame, ToFrame> rotation_;

  template<typename FromFrame, typename ToFrame>
  friend class Permutation;
  template<typename FromFrame, typename ToFrame>
  friend class Rotation;

  template<typename FromFrame, typename ThroughFrame, typename ToFrame>
  friend OrthogonalMap<FromFrame, ToFrame> operator*(
      OrthogonalMap<ThroughFrame, ToFrame> const& left,
      OrthogonalMap<FromFrame, ThroughFrame> const& right);

  friend class OrthogonalMapTest;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "Geometry/orthogonal_map_body.hpp"
