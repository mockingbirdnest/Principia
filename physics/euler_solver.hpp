#pragma once

#include <optional>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/signature.hpp"
#include "geometry/space.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _euler_solver {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_rigid_motion;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A solver for Euler's rotation equations.  It follows [CFSZ07] and [CFSZ08].
// See documentation/Celledoni.pdf for corrections and adaptations.
template<typename InertialFrame, typename PrincipalAxesFrame>
class EulerSolver {
  static_assert(!InertialFrame::may_rotate);

 public:
  using AttitudeRotation = Rotation<PrincipalAxesFrame, InertialFrame>;

  // Constructs a solver for a body with the given moments_of_inertia in its
  // principal axes frame.  The moments must be in increasing order.  At
  // initial_time the angular momentum is initial_angular_momentum and the
  // attitude initial_attitude.
  EulerSolver(
      R3Element<MomentOfInertia> const& moments_of_inertia,
      Bivector<AngularMomentum, InertialFrame> const& initial_angular_momentum,
      AttitudeRotation const& initial_attitude,
      Instant const& initial_time);

  R3Element<MomentOfInertia> const& moments_of_inertia() const;

  // Computes the angular momentum at the given time in the principal axes.
  // This is mostly useful as input to the following two functions.
  Bivector<AngularMomentum, PrincipalAxesFrame> AngularMomentumAt(
      Instant const& time) const;

  AngularVelocity<PrincipalAxesFrame> AngularVelocityFor(
      Bivector<AngularMomentum, PrincipalAxesFrame> const& angular_momentum)
      const;

  // Computes the attitude at the given time, using the angular momentum
  // computed by the previous function.
  AttitudeRotation AttitudeAt(
      Bivector<AngularMomentum, PrincipalAxesFrame> const& angular_momentum,
      Instant const& time) const;

  // Equivalent to this->AttitudeAt(this->AngularMomentumAt(time), time), where
  // the angular momentum is not needed.
  AttitudeRotation AttitudeAt(Instant const& time) const;

  // The motion of the body at the given time.  The centre of gravity of the
  // body moves according to |linear_motion|.
  RigidMotion<PrincipalAxesFrame, InertialFrame> MotionAt(
      Instant const& time,
      DegreesOfFreedom<InertialFrame> const& linear_motion) const;

  void WriteToMessage(not_null<serialization::EulerSolver*> message) const;
  static EulerSolver ReadFromMessage(serialization::EulerSolver const& message);

 private:
  using ℬₜ = Frame<struct ℬₜTag>;
  using ℬʹ = Frame<struct ℬʹTag>;

  // A frame which is rotated from PrincipalAxesFrame such that the coordinates
  // of m along which we project is positive.  Used for all internal
  // computations.
  using PreferredPrincipalAxesFrame =
      Frame<struct PreferredPrincipalAxesFrameTag>;

  using PreferredAngularMomentumBivector =
      Bivector<AngularMomentum, PreferredPrincipalAxesFrame>;

  // The formula to use, following [CFSZ07], Section 2.2.  They don't have a
  // formula for the spherical case.
  enum class Formula {
    i,
    ii,
    iii,
    Sphere,
  };

  // For case (iii) we use project either on e₁ (as we do in case (i)) or on
  // e₃ (as we do in case (ii)) depending on which of the x and z coordinates of
  // m is larger (in absolute value).
  enum class Region {
    e₁,
    e₃,
    Motionless,
  };

  Rotation<PreferredPrincipalAxesFrame, ℬₜ> Compute𝒫ₜ(
      PreferredAngularMomentumBivector const& angular_momentum) const;

  // Construction parameters.
  R3Element<MomentOfInertia> const moments_of_inertia_;
  Bivector<AngularMomentum, InertialFrame> const
      serialized_initial_angular_momentum_;
  AttitudeRotation const initial_attitude_;
  Instant const initial_time_;
  AngularMomentum const G_;
  PreferredAngularMomentumBivector initial_angular_momentum_;
  Rotation<ℬʹ, InertialFrame> ℛ_;

  // A signature that describes which axes are flipped to adjust the signs of
  // the coordinates of m.  It incorporates σ, σʹ and σʺ from [CFSZ07].
  Signature<PrincipalAxesFrame, PreferredPrincipalAxesFrame> 𝒮_;

  // Importantly, the formula and the region to use are constants of motion.
  Formula formula_;
  Region region_;

  // Only the parameters needed for the selected formula are non-NaN after
  // construction.

  AngularFrequency λ_ = NaN<AngularFrequency>;

  AngularMomentum B₂₃_ = NaN<AngularMomentum>;
  AngularMomentum B₁₃_ = NaN<AngularMomentum>;
  AngularMomentum B₃₁_ = NaN<AngularMomentum>;
  AngularMomentum B₂₁_ = NaN<AngularMomentum>;

  double n_ = NaN<double>;
  double mc_ = NaN<double>;
  Angle ν_ = NaN<Angle>;
  Angle ψ_offset_ = NaN<Angle>;
  double ψ_arctan_multiplier_ = NaN<double>;
  MomentOfInertia ψ_cn_multiplier_ = NaN<MomentOfInertia>;
  MomentOfInertia ψ_sn_multiplier_ = NaN<MomentOfInertia>;
  AngularMomentum ψ_cosh_multiplier_ = NaN<AngularMomentum>;
  AngularMomentum ψ_sinh_multiplier_ = NaN<AngularMomentum>;
  double ψ_elliptic_pi_multiplier_ = NaN<double>;
  AngularFrequency ψ_t_multiplier_ = NaN<AngularFrequency>;
};

}  // namespace internal

using internal::EulerSolver;

}  // namespace _euler_solver
}  // namespace physics
}  // namespace principia

#include "physics/euler_solver_body.hpp"
