#pragma once

#include <functional>
#include <type_traits>

#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/conformal_map.hpp"
#include "geometry/frame.hpp"
#include "geometry/homothecy.hpp"
#include "geometry/space_transformations.hpp"
#include "geometry/space.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {

FORWARD_DECLARE_FROM(rigid_motion,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     RigidMotion);

namespace _similar_motion {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_conformal_map;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_homothecy;
using namespace principia::geometry::_space_transformations;
using namespace principia::geometry::_space;
using namespace principia::physics::_rigid_motion;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// The instantaneous motion of |ToFrame| with respect to |FromFrame|.  This is
// the derivative of a |Similarity<FromFrame, ToFrame>|.  In order to invert,
// the |Similarity| is needed, and we need its linear part anyway, so we store
// it (and we forward its action on positions).
template<typename FromFrame, typename ToFrame>
class SimilarMotion final {
  template<typename T>
  using other_frame_t = other_type_t<T, FromFrame, ToFrame>;

 public:
  Similarity<FromFrame, ToFrame> const& similarity() const;
  ConformalMap<double, FromFrame, ToFrame> const& conformal_map() const;

  template<typename F>
  AngularVelocity<other_frame_t<F>> angular_velocity_of() const;
  template<typename F>
  Velocity<other_frame_t<F>> velocity_of_origin_of() const;

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;

  SimilarMotion<ToFrame, FromFrame> Inverse() const;

  static SimilarMotion DilatationAboutOrigin(
      Homothecy<double, FromFrame, ToFrame> const& homothecy,
      Variation<double> dilatation_rate);

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
  SimilarMotion(Similarity<FromFrame, ToFrame> similarity,
                AngularVelocity<FromFrame> angular_velocity_of_to_frame,
                Velocity<FromFrame> velocity_of_to_frame_origin,
                Variation<double> dilatation_rate_of_to_frame);

  Similarity<FromFrame, ToFrame> similarity_;
  AngularVelocity<FromFrame> angular_velocity_of_to_frame_;
  Velocity<FromFrame> velocity_of_to_frame_origin_;
  Variation<double> dilatation_rate_of_to_frame_;

  template<typename From, typename To>
  friend class _rigid_motion::RigidMotion;

  template<typename From, typename To>
  friend class SimilarMotion;

  template<typename From, typename Through, typename To>
  friend SimilarMotion<From, To> operator*(
      SimilarMotion<Through, To> const& left,
      SimilarMotion<From, Through> const& right);

  template<typename, typename>
  friend class RigidMotion;
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
