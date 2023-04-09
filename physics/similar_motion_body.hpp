#pragma once

#include "geometry/identity.hpp"
#include "physics/similar_motion.hpp"

namespace principia {
namespace physics {
namespace _similar_motion {
namespace internal {

using namespace principia::geometry::_rotation;

template<typename FromFrame, typename ToFrame>
inline Similarity<FromFrame, ToFrame> const&
SimilarMotion<FromFrame, ToFrame>::similarity() const {
  return similarity_;
}

template<typename FromFrame, typename ToFrame>
ConformalMap<double, FromFrame, ToFrame> const&
SimilarMotion<FromFrame, ToFrame>::conformal_map() const {
  return similarity_.linear_map();
}

template<typename FromFrame, typename ToFrame>
template<typename F>
AngularVelocity<
    typename SimilarMotion<FromFrame, ToFrame>::template other_frame_t<F>>
SimilarMotion<FromFrame, ToFrame>::angular_velocity_of() const {
  if constexpr (std::is_same_v<F, ToFrame>) {
    return angular_velocity_of_to_frame_;
  } else if constexpr (std::is_same_v<F, FromFrame>) {
    return -similarity_.linear_map()(angular_velocity_of_to_frame_);
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
    return velocity_of_to_frame_origin_;
  } else if constexpr (std::is_same_v<F, FromFrame>) {
    return (*this)({FromFrame::origin, FromFrame::unmoving}).velocity();
  } else {
    static_assert(std::is_same_v<F, ToFrame> || std::is_same_v<F, FromFrame>,
                  "Nonsensical frame");
  }
}

template<typename FromFrame, typename ToFrame>
DegreesOfFreedom<ToFrame> SimilarMotion<FromFrame, ToFrame>::operator()(
    DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const {
  Position<ToFrame> position = similarity_(degrees_of_freedom.position());
  return {position,
          dilatation_rate_of_to_frame_ * (position - ToFrame::origin) +
              angular_velocity_of<FromFrame>() * (position - ToFrame::origin) /
                  Radian +
              conformal_map()(degrees_of_freedom.velocity() -
                              velocity_of_to_frame_origin_)};
}

template<typename FromFrame, typename ToFrame>
SimilarMotion<ToFrame, FromFrame>
SimilarMotion<FromFrame, ToFrame>::Inverse() const {
  return SimilarMotion<ToFrame, FromFrame>(similarity_.Inverse(),
                                           angular_velocity_of<FromFrame>(),
                                           velocity_of_origin_of<FromFrame>(),
                                           -dilatation_rate_of_to_frame_);
}

template<typename FromFrame, typename ToFrame>
SimilarMotion<FromFrame, ToFrame>
SimilarMotion<FromFrame, ToFrame>::DilatationAboutOrigin(
    Homothecy<double, FromFrame, ToFrame> const& homothecy,
    Variation<double> dilatation_rate) {
  return SimilarMotion(
      Similarity<FromFrame, ToFrame>(FromFrame::origin,
                                     ToFrame::origin,
                                     homothecy.template Forget<ConformalMap>()),
      AngularVelocity<FromFrame>{},
      Velocity<FromFrame>{},
      dilatation_rate);
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
SimilarMotion<FromFrame, ToFrame>
SimilarMotion<FromFrame, ToFrame>::Identity() {
  return SimilarMotion(Similarity<FromFrame, ToFrame>::Identity(),
                       AngularVelocity<FromFrame>{},
                       Velocity<FromFrame>{},
                       Variation<double>{});
}

template<typename FromFrame, typename ToFrame>
SimilarMotion<FromFrame, ToFrame>::SimilarMotion(
    Similarity<FromFrame, ToFrame> similarity,
    AngularVelocity<FromFrame> angular_velocity_of_to_frame,
    Velocity<FromFrame> velocity_of_to_frame_origin,
    Variation<double> dilatation_rate_of_to_frame)
    : similarity_(similarity),
      angular_velocity_of_to_frame_(angular_velocity_of_to_frame),
      velocity_of_to_frame_origin_(velocity_of_to_frame_origin),
      dilatation_rate_of_to_frame_(dilatation_rate_of_to_frame) {}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
SimilarMotion<FromFrame, ToFrame> operator*(
    SimilarMotion<ThroughFrame, ToFrame> const& left,
    SimilarMotion<FromFrame, ThroughFrame> const& right) {
  return SimilarMotion<FromFrame, ToFrame>(
      left.similarity() * right.similarity(),
      right.angular_velocity_of_to_frame_ +
          right.conformal_map().Inverse()(left.angular_velocity_of_to_frame_),
      right.Inverse()(left.Inverse()({ToFrame::origin, ToFrame::unmoving}))
          .velocity(),
      left.dilatation_rate_of_to_frame_ + right.dilatation_rate_of_to_frame_);
}

}  // namespace internal
}  // namespace _similar_motion
}  // namespace physics
}  // namespace principia
