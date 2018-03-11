
#pragma once

#include <functional>
#include <vector>

#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/integrators.pb.h"

namespace principia {
namespace integrators {

// The |Solve| function of the |AdaptiveStepSizeIntegrator| exclusively returns
// one of the following statuses.
namespace termination_condition {
constexpr base::Error Done = base::Error::OK;
// The integration may be retried with the same arguments and progress will
// happen.
constexpr base::Error ReachedMaximalStepCount = base::Error::ABORTED;
// A singularity.
constexpr base::Error VanishingStepSize = base::Error::FAILED_PRECONDITION;
}  // namespace termination_condition

namespace internal_ordinary_differential_equations {

using base::Error;
using base::not_null;
using base::Status;
using geometry::Instant;
using numerics::DoublePrecision;
using quantities::Difference;
using quantities::Time;
using quantities::Variation;

// A differential equation of the form y′ = f(y, t).
// This should use parameter packs, but it confuses MSVC...
template<typename... State>
struct FirstOrderOrdinaryDifferentialEquation final {
  using RightHandSideComputation = std::function<Status(
      Instant const& t,
      std::tuple<std::vector<State>...> const& state,
      std::tuple<std::vector<quantities::Quotient<Difference<State, State>, Time>>...>&
          derivatives)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(std::tuple<std::vector<State>...> const& y, Instant const& t);

    std::tuple<std::vector<DoublePrecision<State>>...> y;
    DoublePrecision<Instant> time;

    friend bool operator==(SystemState const& lhs, SystemState const& rhs) {
      return lhs.y == rhs.y && lhs.time == rhs.time;
    }
  };

  using SystemStateError = std::tuple<std::vector<Difference<State, State>>...>;

  // A functor that computes f(y, t) and stores it in |derivatives|.
  // This functor must be called with |std::get<i>(derivatives).size()| equal to
  // |std::get<i>(state).size()| for all i, but there is no requirement on the
  // values in |derivatives|.
  RightHandSideComputation compute_derivative;
};

// TODO(eggrobin): remove the following, this is just testing that it compiles.
FirstOrderOrdinaryDifferentialEquation<quantities::Length,
                                       quantities::Speed,
                                       quantities::Mass,
                                       double> ode;
FirstOrderOrdinaryDifferentialEquation<quantities::Length,
                                       quantities::Speed,
                                       quantities::Mass,
                                       double>::SystemState ode_state;
FirstOrderOrdinaryDifferentialEquation<quantities::Length,
                                       quantities::Speed,
                                       quantities::Mass,
                                       double>::SystemState ode_state_error;

// A differential equation of the form q″ = f(q, t).
// |Position| is the type of q.
template<typename Position_>
struct SpecialSecondOrderDifferentialEquation final {
  using Position = Position_;
  // The type of Δq.
  using Displacement = Difference<Position>;
  // The type of q′.
  using Velocity = Variation<Position>;
  // The type of q″.
  using Acceleration = Variation<Velocity>;
  using RightHandSideComputation =
      std::function<
          Status(Instant const& t,
                 std::vector<Position> const& positions,
                 std::vector<Acceleration>& accelerations)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(std::vector<Position> const& q,
                std::vector<Velocity> const& v,
                Instant const& t);

    std::vector<DoublePrecision<Position>> positions;
    std::vector<DoublePrecision<Velocity>> velocities;
    DoublePrecision<Instant> time;

    friend bool operator==(SystemState const& lhs, SystemState const& rhs) {
      return lhs.positions == rhs.positions &&
             lhs.velocities == rhs.velocities &&
             lhs.time == rhs.time;
    }

    void WriteToMessage(not_null<serialization::SystemState*> message) const;
    static SystemState ReadFromMessage(
        serialization::SystemState const& message);
  };

  struct SystemStateError final {
    std::vector<Displacement> position_error;
    std::vector<Velocity> velocity_error;
  };

  // A functor that computes f(q, t) and stores it in |accelerations|.
  // This functor must be called with |accelerations.size()| equal to
  // |positions.size()|, but there is no requirement on the values in
  // |acceleration|.
  RightHandSideComputation compute_acceleration;
};

// An initial value problem.
template<typename ODE>
struct IntegrationProblem final {
  ODE equation;
  typename ODE::SystemState initial_state;
};

}  // namespace internal_ordinary_differential_equations

using internal_ordinary_differential_equations::IntegrationProblem;
using internal_ordinary_differential_equations::
    SpecialSecondOrderDifferentialEquation;

}  // namespace integrators
}  // namespace principia

#include "integrators/ordinary_differential_equations_body.hpp"
