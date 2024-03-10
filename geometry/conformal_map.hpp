#pragma once

#include "base/concepts.hpp"
#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "geometry/signature.hpp"
#include "geometry/space.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {

FORWARD_DECLARE(TEMPLATE(typename Scalar,
                         typename FromFrame,
                         typename ToFrame) class,
                Homothecy,
                FROM(homothecy),
                INTO(conformal_map));
FORWARD_DECLARE(TEMPLATE(typename FromFrame, typename ToFrame) class,
                OrthogonalMap,
                FROM(orthogonal_map),
                INTO(conformal_map));

namespace _conformal_map {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_mappable;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_linear_map;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;

// This is really a *linear* conformal map, but we don't call it
// LinearConformalMap because that's a mouthful.
template<typename Scalar, typename FromFrame, typename ToFrame>
class ConformalMap : public LinearMap<ConformalMap<Scalar, FromFrame, ToFrame>,
                                      FromFrame, ToFrame> {
 public:
  // The only way to construct conformal maps is as a product of conformal maps
  // obtained by forgetting other linear maps.

  Scalar scale() const;

  Cube<Scalar> Determinant() const;

  ConformalMap<Inverse<Scalar>, ToFrame, FromFrame> Inverse() const;

  template<typename S = Scalar,
           typename = std::enable_if_t<std::is_floating_point_v<S> ||
                                       std::is_integral_v<S>>>
  static ConformalMap Identity();

  template<typename VScalar>
  Vector<Product<VScalar, Scalar>, ToFrame> operator()(
      Vector<VScalar, FromFrame> const& vector) const;

  AngularVelocity<ToFrame> operator()(
      AngularVelocity<FromFrame> const& angular_velocity) const;

  template<typename T>
  typename Mappable<ConformalMap, T>::type operator()(T const& t) const;

  // This is the orthogonal map that applies to multivectors that are once
  // contravariant, once covariant.
  OrthogonalMap<FromFrame, ToFrame> orthogonal_map¹₁() const;

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  static ConformalMap ReadFromMessage(serialization::LinearMap const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

  void WriteToMessage(not_null<serialization::ConformalMap*> message) const;
  static ConformalMap ReadFromMessage(
      serialization::ConformalMap const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

 private:
  ConformalMap(Scalar const& scale,
               Quaternion const& quaternion);

  using RotatedAndSignedFrame = Frame<struct RotatedAndSignedFrameTag,
                                      ToFrame::motion,
                                      ToFrame::handedness>;
  using SignedFrame = Frame<struct SignedFrameTag,
                            ToFrame::motion,
                            ToFrame::handedness>;

  static constexpr Signature<FromFrame, SignedFrame> MakeSignature();
  Rotation<SignedFrame, RotatedAndSignedFrame> MakeRotation() const;
  Homothecy<Scalar, RotatedAndSignedFrame, ToFrame> MakeHomothecy() const;

  Scalar scale_;
  Quaternion quaternion_;

  template<typename S, typename From, typename To>
  friend class ConformalMap;
  template<typename S, typename From, typename To>
  friend class _homothecy::Homothecy;
  template<typename From, typename To>
  friend class _orthogonal_map::OrthogonalMap;

  template<typename L, typename R,
           typename From, typename Through, typename To>
  friend ConformalMap<Product<L, R>, From, To> operator*(
      ConformalMap<L, Through, To> const& left,
      ConformalMap<R, From, Through> const& right);

  template<typename S, typename From, typename To>
  friend std::ostream& operator<<(
      std::ostream& out,
      ConformalMap<S, From, To> const& conformal_map);
};

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
ConformalMap<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    ConformalMap<LScalar, ThroughFrame, ToFrame> const& left,
    ConformalMap<RScalar, FromFrame, ThroughFrame> const& right);

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    ConformalMap<Scalar, FromFrame, ToFrame> const& conformal_map);

}  // namespace internal

using internal::ConformalMap;

}  // namespace _conformal_map
}  // namespace geometry
}  // namespace principia

#include "geometry/conformal_map_body.hpp"
