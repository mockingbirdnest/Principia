#pragma once

#include "geometry/identity.hpp"
#include "physics/similar_motion.hpp"

namespace principia {
namespace physics {
namespace _similar_motion {
namespace internal {

using namespace principia::geometry::_identity;

template<typename FromFrame, typename ToFrame>
using IdentityLinearMap = Identity<FromFrame, ToFrame>;

template<typename FromFrame, typename ToFrame>
template<typename ThroughFrame>
SimilarMotion<FromFrame, ToFrame>::SimilarMotion(
    RigidMotion<FromFrame, ThroughFrame> const& rigid_motion,
    Homothecy<double, ThroughFrame, ToFrame> const& dilatation,
    Variation<double> const& dilatition_rate_of_to_frame)
    : rigid_motion_(RigidMotion<ThroughFrame, Through>::Identity() *
                    rigid_motion),
      dilatation_(dilatation *
                  Homothecy<double, Through, ThroughFrame>::Identity()),
      dilatation_rate_of_to_frame_(dilatition_rate_of_to_frame) {}

template<typename FromFrame, typename ToFrame>
template<typename ThroughFrame>
SimilarMotion<FromFrame, ToFrame>::SimilarMotion(
    Homothecy<double, FromFrame, ThroughFrame> const& dilatation,
    RigidMotion<ThroughFrame, ToFrame> const& rigid_motion,
    Variation<double> const& dilatation_rate_of_to_frame)
    : SimilarMotion(
          RigidMotion<ToFrame, ThroughFrame>::Identity() * rigid_motion *
              RigidMotion<FromFrame, ThroughFrame>::Identity(),
          Homothecy<double, ThroughFrame, ToFrame>::Identity() * dilatation *
              Homothecy<double, ThroughFrame, FromFrame>::Identity(),
          dilatation_rate_of_to_frame) {}

template<typename FromFrame, typename ToFrame>
ConformalMap<double, FromFrame, ToFrame>
SimilarMotion<FromFrame, ToFrame>::conformal_map() const {
  return dilatation_.Forget<ConformalMap>() *
         rigid_motion_.orthogonal_map().Forget<ConformalMap>();
}

template<typename FromFrame, typename ToFrame>
template<typename F>
AngularVelocity<
    typename SimilarMotion<FromFrame, ToFrame>::template other_frame_t<F>>
SimilarMotion<FromFrame, ToFrame>::angular_velocity_of() const {
  if constexpr (std::is_same_v<F, ToFrame>) {
    return rigid_motion_.angular_velocity_of<Through>();
  } else if constexpr (std::is_same_v<F, FromFrame>) {
    return IdentityLinearMap<Through, ToFrame>()(
        rigid_motion_.angular_velocity_of<FromFrame>());
  } else {
    static_assert(std::is_same_v<F, ToFrame> || std::is_same_v<F, FromFrame>,
                  "Nonsensical frame");
  }
}

template<typename FromFrame, typename ToFrame>
template<typename F>
Velocity<typename SimilarMotion<FromFrame, ToFrame>::template other_frame_t<F>>
SimilarMotion<FromFrame, ToFrame>::velocity_of_origin_of() const {
  if constexpr (std::is_same_v<F, ToFrame>) {
    // No use of |dilatation_| on this path because the origin of |Through| is
    // the same as the origin of |ToFrame|.
    return rigid_motion_.velocity_of_origin_of<Through>();
  } else if constexpr (std::is_same_v<F, FromFrame>) {
    return dilatation_(rigid_motion_.velocity_of_origin_of<FromFrame>());
  } else {
    static_assert(std::is_same_v<F, ToFrame> || std::is_same_v<F, FromFrame>,
                  "Nonsensical frame");
  }
}

template<typename FromFrame, typename ToFrame>
DegreesOfFreedom<ToFrame> SimilarMotion<FromFrame, ToFrame>::operator()(
    DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const {
  auto const degrees_of_freedom_in_through = rigid_motion_(degrees_of_freedom);
  auto const& qᴿ = degrees_of_freedom_in_through.position();
  auto const& q̇ᴿ = degrees_of_freedom_in_through.velocity();
  auto const qᴾ = dilatation_(qᴿ - Through::origin);
  return {qᴾ + ToFrame::origin,
          dilatation_rate_of_to_frame_ * qᴾ + dilatation_(q̇ᴿ)};
}

template<typename FromFrame, typename ToFrame>
SimilarMotion<ToFrame, FromFrame>
SimilarMotion<FromFrame, ToFrame>::Inverse() const {
  SimilarMotion<ToFrame, Through> const inverse_dilating_motion(
      RigidMotion<ToFrame, ToFrame>::Identity(),
      dilatation_.Inverse(),
      -dilatation_rate_of_to_frame_);
  return rigid_motion_.Inverse().Forget<SimilarMotion>() *
         inverse_dilating_motion;
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
SimilarMotion<FromFrame, ToFrame>
SimilarMotion<FromFrame, ToFrame>::Identity() {
  return SimilarMotion(RigidMotion<FromFrame, Through>::Identity(),
                       Homothecy<double, Through, ToFrame>::Identity(),
                       Variation<double>{});
}

template<typename FromFrame, typename ToFrame>
template<typename ThroughOtherFrame>
SimilarMotion<FromFrame, ToFrame>::CommutedSplit<
    ThroughOtherFrame>::CommutedSplit(SimilarMotion const& similar_motion)
    : dilatation(Homothecy<double, Through, ThroughOtherFrame>::Identity() *
                 similar_motion.dilatation_.Inverse()),
      rigid_motion(similar_motion.rigid_motion_.Inverse() *
                   RigidMotion<ThroughOtherFrame, Through>::Identity()) {}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
SimilarMotion<FromFrame, ToFrame> operator*(
    SimilarMotion<ThroughFrame, ToFrame> const& left,
    SimilarMotion<FromFrame, ThroughFrame> const& right) {
  using Left = SimilarMotion<ThroughFrame, ToFrame>;
  using Right = SimilarMotion<FromFrame, ThroughFrame>;
  typename Left::template CommutedSplit<typename Right::Through> split(left);
  auto const& hL = split.dilatation;
  auto const& rL = split.rigid_motion;
  auto const& rR = right.rigid_motion_;
  auto const& hR = right.dilatation_;
  auto const r = rL * rR;
}

}  // namespace internal
}  // namespace _similar_motion
}  // namespace physics
}  // namespace principia
