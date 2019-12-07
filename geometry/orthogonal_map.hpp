
#pragma once

#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

namespace physics {
class RigidMotionTest;
}  // namespace physics

namespace geometry {

FORWARD_DECLARE_FROM(identity,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     Identity);
FORWARD_DECLARE_FROM(permutation,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     Permutation);
FORWARD_DECLARE_FROM(rotation,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     Rotation);
FORWARD_DECLARE_FROM(symmetric_bilinear_form,
                     TEMPLATE(typename Scalar, typename Frame) class,
                     SymmetricBilinearForm);

namespace internal_orthogonal_map {

using base::not_null;

// An orthogonal map between the inner product spaces |FromFrame| and
// |ToFrame|, as well as the induced maps on the exterior algebra.
// The orthogonal map is modeled as a rotoinversion.
template<typename FromFrame, typename ToFrame>
class OrthogonalMap : public LinearMap<FromFrame, ToFrame> {
 public:
  Sign Determinant() const override;

  Rotation<FromFrame, ToFrame> const& rotation() const;

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

  template<typename Scalar>
  SymmetricBilinearForm<Scalar, ToFrame> operator()(
      SymmetricBilinearForm<Scalar, FromFrame> const& form) const;

  template<typename T>
  typename base::Mappable<OrthogonalMap, T>::type operator()(T const& t) const;

  static OrthogonalMap Identity();

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  static OrthogonalMap ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::OrthogonalMap*> message) const;
  static OrthogonalMap ReadFromMessage(
      serialization::OrthogonalMap const& message);

 private:
  OrthogonalMap(Sign const& determinant,
                Rotation<FromFrame, ToFrame> const& rotation);

  Sign determinant_;
  Rotation<FromFrame, ToFrame> rotation_;

  template<typename From, typename To>
  friend class internal_identity::Identity;
  template<typename From, typename To>
  friend class OrthogonalMap;
  template<typename From, typename To>
  friend class internal_permutation::Permutation;
  template<typename From, typename To>
  friend class internal_rotation::Rotation;

  template<typename From, typename Through, typename To>
  friend OrthogonalMap<From, To> operator*(
      OrthogonalMap<Through, To> const& left,
      OrthogonalMap<From, Through> const& right);

  friend class OrthogonalMapTest;

  // TODO(phl): This friendship could be avoided if we had symmetries.
  friend class physics::RigidMotionTest;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map);

}  // namespace internal_orthogonal_map

using internal_orthogonal_map::OrthogonalMap;

}  // namespace geometry
}  // namespace principia

#include "geometry/orthogonal_map_body.hpp"
