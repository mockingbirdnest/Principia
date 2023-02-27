#pragma once

#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/plane.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/integration_parameters.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using integrators::AdaptiveStepSizeIntegrator;
using integrators::ExplicitFirstOrderOrdinaryDifferentialEquation;
using quantities::Acceleration;
using quantities::Angle;
using quantities::Difference;
using quantities::Infinity;
using quantities::Length;
using quantities::SpecificEnergy;
using quantities::si::Metre;
using quantities::si::Second;
using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_named_quantities;
using namespace principia::geometry::_plane;

template<typename InertialFrame, typename Frame>
class Equipotential {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  // Equipotential lines are parameterized by this dimentionless quantity.
  using IndependentVariable = double;

  // The first state variable is a point of an equipotential.  The second state
  // variable is the (scale-free) intensity of the braking.
  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<IndependentVariable,
                                                     Position<Frame>,
                                                     double>;
  using DependentVariables = typename ODE::DependentVariables;
  using AdaptiveParameters = physics::AdaptiveStepParameters<ODE>;
  using Line = std::vector<DependentVariables>;
  using Lines = std::vector<Line>;

  Equipotential(
      AdaptiveParameters const& adaptive_parameters,
      not_null<DynamicFrame<InertialFrame, Frame> const*> dynamic_frame);

  // Computes an equipotential line going through the given point.
  Line ComputeLine(Plane<Frame> const& plane,
                   Instant const& t,
                   Position<Frame> const& position) const;

  // Computes an equipotential line for the total energy determined by the
  // |degrees_of_freedom|.
  Line ComputeLine(Plane<Frame> const& plane,
                   Instant const& t,
                   DegreesOfFreedom<Frame> const& degrees_of_freedom) const;

  // Computes an equipotential line for the given |total_energy| starting from
  // |start_position|.
  Line ComputeLine(Plane<Frame> const& plane,
                   Instant const& t,
                   Position<Frame> const& start_position,
                   SpecificEnergy const& total_energy) const;

  // Computes equipotential lines for the given |total_energy|.  Each of the
  // given |start_positions| ends up enclosed by exactly one line of the result.
  // The |start_positions| must be coplanar in a plane parallel to |plane|.
  Lines ComputeLines(Plane<Frame> const& plane,
                     Instant const& t,
                     std::vector<Position<Frame>> const& start_positions,
                     SpecificEnergy const& total_energy) const;

  struct Well {
    Position<Frame> position;
    // Below this radius the potential has spurious maxima (and eventually NaNs
    // and infinities), so algorithms should not look there.
    Length radius;
  };

  // Computes equipotential lines for the given |energy| that delineate the
  // |peaks| from the |wells| and from the “well at infinity” (which is a well
  // because we are in a rotating frame).  An equipotential delineates a peak
  // from a well if it encloses the peak but not the well, or vice-versa.  Given
  // a position, |towards_infinity| should return a position far away where the
  // potential is lower, in a direction where not much happens, e.g., away from
  // the centre in a rotating frame.
  Lines ComputeLines(
      Plane<Frame> const& plane,
      Instant const& t,
      std::vector<Position<Frame>> const& peaks,
      std::vector<Well> const& wells,
      std::function<Position<Frame>(Position<Frame>)> towards_infinity,
      SpecificEnergy const& energy) const;

 private:
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;
  using DependentVariableDerivatives =
      typename ODE::DependentVariableDerivatives;
  using State = typename ODE::State;

  static constexpr IndependentVariable const s_initial_ = 0;
  static constexpr IndependentVariable const s_final_ =
      Infinity<IndependentVariable>;
  static constexpr IndependentVariableDifference const initial_s_step_ = 1;
  static constexpr Length const characteristic_length_ = 1 * Metre;

  // TODO(phl): One or both of these values should probably be a parameter.
  static constexpr double β_max_ = 1e6;
  static constexpr double β_tolerance_ = 1;

  // The |binormal| determines in what direction we go around the curve.  It may
  // be anything, but must be consistent across calls to the right-hand side.
  absl::Status RightHandSide(Bivector<double, Frame> const& binormal,
                             Position<Frame> const& position,
                             Instant const& t,
                             IndependentVariable s,
                             DependentVariables const& values,
                             DependentVariableDerivatives& derivatives) const;

  double ToleranceToErrorRatio(IndependentVariableDifference current_s_step,
                               State const& /*state*/,
                               typename State::Error const& error) const;

  // Computes the winding number of |line| around |position|.  |line| and
  // |position| must be in a plane paralles to |plane|.  The returned integer is
  // nonnegative, i.e., doesn't give information about the direction in which
  // the |line| rotates around |position|.
  std::int64_t WindingNumber(Plane<Frame> const& plane,
                             Position<Frame> const& position,
                             std::vector<Position<Frame>> const& line) const;

  AdaptiveParameters const& adaptive_parameters_;
  not_null<DynamicFrame<InertialFrame, Frame> const*> const dynamic_frame_;
};

}  // namespace internal_equipotential

using internal_equipotential::Equipotential;

}  // namespace physics
}  // namespace principia

#include "physics/equipotential_body.hpp"
