#pragma once

#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

// An orthogonal map between the inner product spaces |FromFrame| and
// |ToFrame|, as well as the induced maps on the exterior algebra.
// The orthogonal map is modeled as a rotoinversion.
template<typename FromFrame, typename ToFrame>
class OrthogonalMap : public LinearMap<FromFrame, ToFrame> {
 public:
  ~OrthogonalMap() override = default;

  Sign Determinant() const override;

  Rotation<FromFrame, ToFrame> rotation() const;

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

  template<typename T>
  typename base::Mappable<OrthogonalMap, T>::type operator()(T const& t) const;

  static OrthogonalMap Identity();

  void WriteToMessage(not_null<serialization::LinearMap*> const message) const;
  static OrthogonalMap ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(
      not_null<serialization::OrthogonalMap*> const message) const;
  static OrthogonalMap ReadFromMessage(
      serialization::OrthogonalMap const& message);

 private:
  OrthogonalMap(Sign const& determinant,
                Rotation<FromFrame, ToFrame> const& rotation);

  Sign determinant_;
  Rotation<FromFrame, ToFrame> rotation_;

  template<typename From, typename To>
  friend class Identity;
  template<typename From, typename To>
  friend class OrthogonalMap;
  template<typename From, typename To>
  friend class Permutation;
  template<typename From, typename To>
  friend class Rotation;

  template<typename From, typename Through, typename To>
  friend OrthogonalMap<From, To> operator*(
      OrthogonalMap<Through, To> const& left,
      OrthogonalMap<From, Through> const& right);

  friend class OrthogonalMapTest;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/orthogonal_map_body.hpp"
