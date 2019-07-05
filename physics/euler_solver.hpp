
#pragma once

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace physics {
namespace internal_euler_solver {

using geometry::Bivector;
using geometry::Frame;
using geometry::Instant;
using quantities::AngularMomentum;
using quantities::Energy;
using quantities::MomentOfInertia;

class EulerSolver {
 public:
  using PrincipalAxesFrame = Frame<serialization::Frame::PhysicsTag,
                                   serialization::Frame::PRINCIPAL_AXES,
                                   /*frame_is_inertial*/ false>;

  EulerSolver(MomentOfInertia const& moment_of_inertia₁,
              MomentOfInertia const& moment_of_inertia₂,
              MomentOfInertia const& moment_of_inertia₃,
              Energy const& kinetic_energy);

  Bivector<AngularMomentum, PrincipalAxesFrame> ComputeAngularMomentum(
      Instant const& t);

};

}  // namespace internal_euler_solver

using internal_euler_solver::EulerSolver;

}  // namespace physics
}  // namespace principia
