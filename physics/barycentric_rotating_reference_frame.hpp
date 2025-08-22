// The files containing the tree of child classes of `RigidReferenceFrame` must
// be included in the order of inheritance to avoid circular dependencies.  This
// class will end up being reincluded as part of the implementation of its
// parent.
#ifndef PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
#include "physics/rigid_reference_frame.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_REFERENCE_FRAME_HPP_
#define PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_REFERENCE_FRAME_HPP_

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/arithmetic.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace physics {
namespace _barycentric_rotating_reference_frame {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_rigid_motion;
using namespace principia::physics::_rigid_reference_frame;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_tuples;

// The origin of the frame is the barycentre of the system.  The X axis
// points to the barycentre of the secondaries.  The Y axis is in the direction
// of the velocity of the barycentre of the secondaries with respect to the
// primary.  The Z axis is in the direction of the angular momentum of the
// system.  The basis has the same orientation as `InertialFrame`.
// Note that if the angular momentum of the system varies, the angular velocity
// of the frame need not be along its Z axis.
template<typename InertialFrame, typename ThisFrame>
class BarycentricRotatingReferenceFrame
    : public RigidReferenceFrame<InertialFrame, ThisFrame> {
  static_assert(ThisFrame::may_rotate);

 public:
  BarycentricRotatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      not_null<MassiveBody const*> primary,
      not_null<MassiveBody const*> secondary);

  BarycentricRotatingReferenceFrame(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      std::vector<not_null<MassiveBody const*>> primaries,
      std::vector<not_null<MassiveBody const*>> secondaries);

  std::vector<not_null<MassiveBody const*>> const& primaries() const;
  std::vector<not_null<MassiveBody const*>> const& secondaries() const;

  template<int degree>
  Derivative<Position<InertialFrame>, Instant, degree> PrimaryDerivative(
      Instant const& t) const;
  template<int degree>
  Derivative<Position<InertialFrame>, Instant, degree> SecondaryDerivative(
      Instant const& t) const;

  Instant t_min() const override;
  Instant t_max() const override;

  RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const override;

  void WriteToMessage(
      not_null<serialization::ReferenceFrame*> message) const override;

  static not_null<std::unique_ptr<BarycentricRotatingReferenceFrame>>
  ReadFromMessage(
      not_null<Ephemeris<InertialFrame> const*> ephemeris,
      serialization::BarycentricRotatingReferenceFrame const& message);

 private:
  using Base = RigidReferenceFrame<InertialFrame, ThisFrame>;
  using BodiesToDegreesOfFreedom =
      typename Ephemeris<InertialFrame>::BodiesToDegreesOfFreedom;
  using BodiesToPositions =
      typename Ephemeris<InertialFrame>::BodiesToPositions;

  template<typename SF, typename SB, int o = 0>
  using Trihedron = typename Base::template Trihedron<SF, SB, o>;

  struct CachedDerivatives {
    Derivatives<Position<InertialFrame>, Instant, 4> derivatives;
    std::array<Instant, 4> times = {Instant() + NaN<Time>,
                                    Instant() + NaN<Time>,
                                    Instant() + NaN<Time>,
                                    Instant() + NaN<Time>};
  };

  template<int degree>
  Derivative<Position<InertialFrame>, Instant, degree> PrimaryDerivative(
      BodiesToDegreesOfFreedom const* bodies_to_degrees_of_freedom,
      Instant const& t) const;
  template<int degree>
  Derivative<Position<InertialFrame>, Instant, degree> PrimaryDerivative(
      BodiesToPositions const* bodies_to_positions,
      Instant const& t) const;
  template<int degree>
  Derivative<Position<InertialFrame>, Instant, degree> SecondaryDerivative(
      BodiesToDegreesOfFreedom const* bodies_to_degrees_of_freedom,
      Instant const& t) const;
  template<int degree>
  Derivative<Position<InertialFrame>, Instant, degree> SecondaryDerivative(
      BodiesToPositions const* bodies_to_positions,
      Instant const& t) const;
  template<
      int degree,
      std::vector<not_null<MassiveBody const*>> const
          BarycentricRotatingReferenceFrame<InertialFrame, ThisFrame>::*bodies>
  Derivative<Position<InertialFrame>, Instant, degree> BarycentreDerivative(
      BodiesToDegreesOfFreedom const* bodies_to_degrees_of_freedom,
      BodiesToPositions const* bodies_to_positions,
      Instant const& t,
      CachedDerivatives& cache) const;

  Vector<Acceleration, InertialFrame> GravitationalAcceleration(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  SpecificEnergy GravitationalPotential(
      Instant const& t,
      Position<InertialFrame> const& q) const override;
  AcceleratedRigidMotion<InertialFrame, ThisFrame> MotionOfThisFrame(
      Instant const& t) const override;

  // Implementation helper that avoids evaluating the degrees of freedom and the
  // accelerations multiple times.
  RigidMotion<InertialFrame, ThisFrame> ToThisFrame(
      Derivatives<Position<InertialFrame>, Instant, 3> const&
          primary_derivative,
      Derivatives<Position<InertialFrame>, Instant, 3> const&
          secondary_derivative) const;

  not_null<Ephemeris<InertialFrame> const*> const ephemeris_;
  std::vector<not_null<MassiveBody const*>> const primaries_;
  std::vector<not_null<MassiveBody const*>> const secondaries_;
  GravitationalParameter const primary_gravitational_parameter_;
  GravitationalParameter const secondary_gravitational_parameter_;
  mutable absl::Mutex lock_;
  // These members optimize costly computations from `BarycentreDerivative` in
  // the frequent case where properties of the frame are repeatedly requested
  // for the same time.
  mutable CachedDerivatives last_evaluated_primary_derivatives_
      GUARDED_BY(lock_);
  mutable CachedDerivatives last_evaluated_secondary_derivatives_
      GUARDED_BY(lock_);
};

}  // namespace internal

using internal::BarycentricRotatingReferenceFrame;

}  // namespace _barycentric_rotating_reference_frame
}  // namespace physics
}  // namespace principia

#include "physics/barycentric_rotating_reference_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BARYCENTRIC_ROTATING_REFERENCE_FRAME_HPP_
#endif  // PRINCIPIA_PHYSICS_RIGID_REFERENCE_FRAME_HPP_
