#pragma once

#include "base/mappable.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
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
FORWARD_DECLARE_FROM(signature,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     Signature);
FORWARD_DECLARE_FROM(
    symmetric_bilinear_form,
    TEMPLATE(typename Scalar,
            typename Frame,
            template<typename, typename> typename Multivector) class,
    SymmetricBilinearForm);

class ConformalMapTest;
class OrthogonalMapTest;

namespace _orthogonal_map {
namespace internal {

using namespace principia::base::_mappable;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_permutation;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_symmetric_bilinear_form;

// An orthogonal map between the inner product spaces |FromFrame| and
// |ToFrame|, as well as the induced maps on the exterior algebra.
// The orthogonal map is modeled as a rotoinversion.
template<typename FromFrame, typename ToFrame>
class OrthogonalMap : public LinearMap<OrthogonalMap<FromFrame, ToFrame>,
                                       FromFrame, ToFrame> {
 public:
  Sign Determinant() const;

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<F::handedness == T::handedness>>
  Rotation<FromFrame, ToFrame> AsRotation() const;

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

  template<typename Scalar,
           template<typename, typename> typename Multivector>
  SymmetricBilinearForm<Scalar, ToFrame, Multivector> operator()(
      SymmetricBilinearForm<Scalar, FromFrame, Multivector> const& form) const;

  template<typename T>
  typename Mappable<OrthogonalMap, T>::type operator()(T const& t) const;

  template<template<typename, typename> typename ConformalMap>
  ConformalMap<FromFrame, ToFrame> Forget() const;
  template<template<typename, typename, typename> typename ConformalMap>
  ConformalMap<double, FromFrame, ToFrame> Forget() const;

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<F::handedness == T::handedness>>
  static OrthogonalMap Identity();

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<is_serializable_v<F> &&
                                       is_serializable_v<T>>>
  static OrthogonalMap ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::OrthogonalMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<is_serializable_v<F> &&
                                       is_serializable_v<T>>>
  static OrthogonalMap ReadFromMessage(
      serialization::OrthogonalMap const& message);

 private:
  explicit OrthogonalMap(Quaternion const& quaternion);

  using IntermediateFrame = Frame<struct IntermediateFrameTag,
                                  ToFrame::motion,
                                  ToFrame::handedness>;

  static constexpr Signature<FromFrame, IntermediateFrame> MakeSignature();
  Rotation<IntermediateFrame, ToFrame> MakeRotation() const;

  Quaternion quaternion_;

  static constexpr Sign determinant_ =
      FromFrame::handedness == ToFrame::handedness ? Sign::Positive()
                                                   : Sign::Negative();

  template<typename From, typename To>
  friend class OrthogonalMap;
  template<typename From, typename To>
  friend class _identity::Identity;
  template<typename From, typename To>
  friend class _permutation::Permutation;
  template<typename From, typename To>
  friend class _rotation::Rotation;
  template<typename From, typename To>
  friend class _signature::Signature;

  template<typename From, typename Through, typename To>
  friend OrthogonalMap<From, To> operator*(
      OrthogonalMap<Through, To> const& left,
      OrthogonalMap<From, Through> const& right);

  template<typename From, typename To>
  friend std::ostream& operator<<(
      std::ostream& out,
      OrthogonalMap<From, To> const& orthogonal_map);

  friend class geometry::ConformalMapTest;
  friend class geometry::OrthogonalMapTest;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map);

}  // namespace internal

using internal::OrthogonalMap;

}  // namespace _orthogonal_map
}  // namespace geometry
}  // namespace principia

namespace principia::geometry {
using namespace principia::geometry::_orthogonal_map;
}  // namespace principia::geometry

#include "geometry/orthogonal_map_body.hpp"
