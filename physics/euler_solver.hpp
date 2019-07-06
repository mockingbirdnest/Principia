
#pragma once

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using geometry::Bivector;
using geometry::Frame;
using geometry::Instant;
using geometry::R3Element;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::AngularMomentum;
using quantities::MomentOfInertia;
using quantities::NaN;

class EulerSolver {
 public:
  using PrincipalAxesFrame = Frame<serialization::Frame::PhysicsTag,
                                   serialization::Frame::PRINCIPAL_AXES,
                                   /*frame_is_inertial*/ false>;
  using AngularMomentumBivector = Bivector<AngularMomentum, PrincipalAxesFrame>;

  EulerSolver(R3Element<MomentOfInertia> const& moments_of_inertia,
              AngularMomentumBivector const& initial_angular_momentum,
              Instant const& initial_time);

  AngularMomentumBivector AngularMomentumAt(Instant const& time) const;

private:
  enum class Formula {
    i,
    ii,
    iii
  };

  Instant const initial_time_;

  // Amusingly, the formula to use is a constant of motion.
  Formula formula_;

  AngularMomentum B₁₃_ = NaN<AngularMomentum>();
  AngularMomentum B₃₁_ = NaN<AngularMomentum>();
  AngularMomentum B₂₁_ = NaN<AngularMomentum>();
  AngularMomentum B₂₃_ = NaN<AngularMomentum>();
  AngularMomentum G_ = NaN<AngularMomentum>();
  AngularFrequency λ₁_ = NaN<AngularFrequency>();
  AngularFrequency λ₃_ = NaN<AngularFrequency>();
  double mc_ = NaN<double>();
  Angle ν_ = NaN<Angle>();
};

}  // namespace internal_euler_solver

using internal_euler_solver::EulerSolver;

}  // namespace physics
}  // namespace principia
