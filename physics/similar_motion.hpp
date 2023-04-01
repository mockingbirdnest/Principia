#pragma once

#include <functional>
#include <type_traits>

#include "base/not_null.hpp"
#include "geometry/affine_map.hpp"
#include "geometry/conformal_map.hpp"
#include "geometry/rigid_transformation.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _similar_motion {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_affine_map;
using namespace principia::geometry::_conformal_map;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_rigid_transformation;
using namespace principia::geometry::_space;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// A trait to determine if Frame is FromFrame or ToFrame and return the other
// one as the |type| member.
//TODO(phl)move to a common place
template<typename Frame, typename FromFrame, typename ToFrame>
struct other_frame;

template<typename Frame>
struct other_frame<Frame, Frame, Frame> {
  using type = Frame;
};

template<typename FromFrame, typename ToFrame>
struct other_frame<FromFrame, FromFrame, ToFrame> {
  using type = ToFrame;
};

template<typename FromFrame, typename ToFrame>
struct other_frame<ToFrame, FromFrame, ToFrame> {
  using type = FromFrame;
};

// The instantaneous motion of |ToFrame| with respect to |FromFrame|.  This is
// the derivative of a |Similarity<FromFrame, ToFrame>|.  In order to invert,
// the |Similarity| is needed, and we need its linear part anyway, so we store
// it (and we forward its action on positions).
template<typename FromFrame, typename ToFrame>
class SimilarMotion final {
  template<typename T>
  using other_frame_t = typename other_frame<T, FromFrame, ToFrame>::type;

 public:
  SimilarMotion(
      Similarity<FromFrame, ToFrame> const& similarity,
      AngularVelocity<FromFrame> const& angular_velocity_of_to_frame,
      Velocity<FromFrame> const& velocity_of_to_frame_origin,
      Variation<double> const& dilatation_rate);

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same_v<F, T>>>
  SimilarMotion(
      Similarity<FromFrame, ToFrame> const& similarity,
      AngularVelocity<ToFrame> const& angular_velocity_of_from_frame,
      Velocity<ToFrame> const& velocity_of_from_frame_origin,
      Variation<double> const& dilatation_rate);

  Similarity<FromFrame, ToFrame> const& similarity() const;
  // Returns |similarity().linear_map()|.
  ConformalMap<double, FromFrame, ToFrame> const& conformal_map() const;

  template<typename F>
  AngularVelocity<other_frame_t<F>> angular_velocity_of() const;
  template<typename F>
  Velocity<other_frame_t<F>> velocity_of_origin_of() const;

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;

  SimilarMotion<ToFrame, FromFrame> Inverse() const;

  // A similar motion expressing that |FromFrame| and |ToFrame| have the same
  // axes, origin, and instantaneous motion.
  // This function is enabled only if both frames have the same handedness (this
  // is a requirement of OrthogonalMap::Identity) and if the |motion| of
  // FromFrame is a special case of that of |ToFrame| (see the comments on
  // |FrameMotion|).
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<(F::handedness == T::handedness &&
                                        F::motion <= T::motion)>>
  static SimilarMotion Identity();

 private:
  Similarity<FromFrame, ToFrame> similarity_;
  // d/dt similarity⁻¹(basis of ToFrame). The positively oriented
  // orthogonal bases of |FromFrame| are acted upon faithfully and transitively
  // by SO(FromFrame), so this lies in the tangent space, i.e., the Lie algebra
  // 𝖘𝔬(FromFrame) ≅ FromFrame ∧ FromFrame.
  AngularVelocity<FromFrame> angular_velocity_of_to_frame_;
  // d/dt similarity⁻¹(ToFrame::origin).
  Velocity<FromFrame> velocity_of_to_frame_origin_;
  Variation<double> dilatation_rate_;

  template<typename, typename>
  friend class SimilarMotion;

  template<typename From, typename Through, typename To>
  friend SimilarMotion<From, To> operator*(
      SimilarMotion<Through, To> const& left,
      SimilarMotion<From, Through> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
SimilarMotion<FromFrame, ToFrame> operator*(
    SimilarMotion<ThroughFrame, ToFrame> const& left,
    SimilarMotion<FromFrame, ThroughFrame> const& right);

}  // namespace internal

using internal::SimilarMotion;

}  // namespace _similar_motion
}  // namespace physics
}  // namespace principia

#include "physics/similar_motion_body.hpp"
