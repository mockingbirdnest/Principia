#pragma once

#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "physics/harmonic_damping.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _geopotential {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_polynomial;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// Representation of the geopotential model of an oblate body.
template<typename Frame>
class Geopotential {
 public:
  // Spherical harmonics will not be damped if their contribution to the radial
  // force exceeds |tolerance| times the central force.
  Geopotential(not_null<OblateBody<Frame> const*> body,
               double tolerance);

  Vector<Quotient<Acceleration, GravitationalParameter>, Frame>
  GeneralSphericalHarmonicsAcceleration(
      Instant const& t,
      Displacement<Frame> const& r,
      Length const& r_norm,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

  Quotient<SpecificEnergy, GravitationalParameter>
  GeneralSphericalHarmonicsPotential(
      Instant const& t,
      Displacement<Frame> const& r,
      Length const& r_norm,
      Square<Length> const& r²,
      Exponentiation<Length, -3> const& one_over_r³) const;

  std::vector<HarmonicDamping> const& degree_damping() const;
  HarmonicDamping const& sectoral_damping() const;

 private:
  // The frame of the surface of the celestial.
  using SurfaceFrame = geometry::_frame::Frame<struct SurfaceFrameTag>;
  static const Vector<double, SurfaceFrame> x_;
  static const Vector<double, SurfaceFrame> y_;

  // These are the types that we return, so better have a name for them.
  using ReducedAcceleration = Quotient<Acceleration, GravitationalParameter>;
  using ReducedPotential = Quotient<SpecificEnergy, GravitationalParameter>;

  // List of reduced quantities computed for all degrees or orders.
  template<int size>
  using ReducedAccelerations =
      std::array<Vector<ReducedAcceleration, Frame>, size>;
  template<int size>
  using ReducedPotentials = std::array<ReducedPotential, size>;

  using UnitVector = Vector<double, Frame>;

  // Holds precomputed data for one evaluation of the acceleration.
  struct Precomputations;

  // Helper templates for iterating over the degrees/orders of the geopotential.
  template<int degree, int order>
  class DegreeNOrderM;
  template<int degree, typename>
  class DegreeNAllOrders;
  template<typename>
  class AllDegrees;

  // |limiting_degree| is the first degree such that
  // |r_norm >= degree_damping_[limiting_degree].outer_threshold()|, or is
  // |degree_damping_.size()| if |r_norm| is below all thresholds.
  // Since |degree_damping_[0].outer_threshold()| and
  // |degree_damping_[1].outer_threshold()| are infinite, |limiting_degree > 1|.
  int LimitingDegree(Length const& r_norm) const;

  not_null<OblateBody<Frame> const*> body_;

  // The contribution from the harmonics of degree n is damped by
  // degree_damping_[n].
  // degree_damping_[0] and degree_damping_[1] have infinite thresholds, and are
  // not used (this class does not compute the central force and disregards
  // degree 1, which is equivalent to a translation of the centre of mass).
  std::vector<HarmonicDamping> degree_damping_;

  // The contribution of the degree 2 sectoral harmonics is damped by
  // |sectoral_damping_|; |degree_damping_[2]| affects only J2.
  // The monotonicity relation
  //   degree_damping[2] ≼ sectoral_damping_ ≼ degree_damping[3]
  // holds, where ≼ denotes the ordering of the thresholds.
  HarmonicDamping sectoral_damping_;
};

}  // namespace internal

using internal::Geopotential;

}  // namespace _geopotential
}  // namespace physics
}  // namespace principia

#include "physics/geopotential_body.hpp"
