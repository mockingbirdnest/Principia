// The files containing the tree of child classes of `ReferenceFrame` must be
// included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
#include "physics/reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_ROTATING_PULSATING_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_ROTATING_PULSATING_REFERENCE_FRAME_HPP_

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "physics/barycentric_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/similar_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/tuples.hpp"

#include "physics/body_centred_body_direction_reference_frame.hpp"

namespace principia {
namespace physics {
namespace _rotating_pulsating_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::physics::_barycentric_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_similar_motion;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_tuples;

template<typename InertialFrame, typename ThisFrame>
class RotatingPulsatingReferenceFrame
    : public ReferenceFrame<InertialFrame, ThisFrame> {
  static_assert(ThisFrame::may_rotate);
  // TODO(egg): This frame may dilate, too.

 public:
  RotatingPulsatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> primary,
      not_null<MassiveBody const*> secondary);

  RotatingPulsatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      std::vector<not_null<MassiveBody const*>> primaries,
      std::vector<not_null<MassiveBody const*>> secondaries);

  std::vector<not_null<MassiveBody const*>> const& primaries() const;
  std::vector<not_null<MassiveBody const*>> const& secondaries() const;

  Instant t_min() const override;
  Instant t_max() const override;

  SimilarMotion<InertialFrame, ThisFrame> ToThisFrameAtTimeSimilarly(
      Instant const& t) const override;

  Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const override;

  Vector<Acceleration, ThisFrame> RotationFreeGeometricAccelerationAtRest(
      Instant const& t,
      Position<ThisFrame> const& position) const override;

  SpecificEnergy GeometricPotential(
      Instant const& t,
      Position<ThisFrame> const& position) const override;

  void WriteToMessage(
      not_null<serialization::ReferenceFrame*> message) const override;

  static not_null<std::unique_ptr<RotatingPulsatingReferenceFrame>>
  ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::RotatingPulsatingReferenceFrame const& message);

 private:
  using RotatingFrame =
      Frame<struct RotatingTag, Arbitrary, ThisFrame::handedness>;

  template<int degree>
  Derivatives<Length, Instant, degree + 1> r_derivatives(
      Instant const& t) const;

  SimilarMotion<ThisFrame, RotatingFrame> ToRotatingFrame(
      Derivatives<Length, Instant, 2> const& r_derivatives_1) const;

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  std::vector<not_null<MassiveBody const*>> const primaries_;
  std::vector<not_null<MassiveBody const*>> const secondaries_;
  BarycentricRotatingReferenceFrame<InertialFrame, RotatingFrame> const
      rotating_frame_;
};

}  // namespace internal

using internal::RotatingPulsatingReferenceFrame;

}  // namespace _rotating_pulsating_reference_frame
}  // namespace physics
}  // namespace principia

#include "physics/rotating_pulsating_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_ROTATING_PULSATING_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_REFERENCE_FRAME_HPP_
