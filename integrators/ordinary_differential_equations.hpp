#pragma once

#include <functional>
#include <tuple>
#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
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

constexpr absl::StatusCode Done = absl::StatusCode::kOk;
// The integration may be retried with the same arguments and progress will
// happen.
constexpr absl::StatusCode ReachedMaximalStepCount = absl::StatusCode::kAborted;
// A singularity.
constexpr absl::StatusCode VanishingStepSize =
    absl::StatusCode::kFailedPrecondition;

// Same as absl::Status::Update, but prefers kAbort.
void UpdateWithAbort(absl::Status const& updater, absl::Status& updated);

}  // namespace termination_condition

namespace internal_ordinary_differential_equations {

using base::not_null;
using geometry::Instant;
using numerics::DoublePrecision;
using quantities::Difference;
using quantities::Quotient;
using quantities::Time;
using quantities::Variation;

// A differential equation of the form y′ = f(s, y).
// |State| is the type of y.
template<typename IndependentVariable_, typename... StateElements>
struct ExplicitFirstOrderOrdinaryDifferentialEquation final {
  static constexpr std::int64_t order = 1;
  using IndependentVariable = IndependentVariable_;
  using IndependentVariableDifference = Difference<IndependentVariable>;
  using State = std::tuple<std::vector<StateElements>...>;
  using StateDifference = std::tuple<std::vector<Difference<StateElements>>...>;
  using StateVariation = std::tuple<std::vector<
      Quotient<Difference<StateElements>, IndependentVariableDifference>>...>;

  using RightHandSideComputation =
      std::function<absl::Status(IndependentVariable const& s,
                                 State const& state,
                                 StateVariation& variations)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(IndependentVariable const& s, State const& y);

    DoublePrecision<IndependentVariable> s;
    std::tuple<std::vector<DoublePrecision<StateElements>>...> y;

    friend bool operator==(SystemState const& lhs, SystemState const& rhs) {
      return lhs.y == rhs.y && lhs.s == rhs.s;
    }

    void WriteToMessage(not_null<serialization::SystemState*> message) const;
    static SystemState ReadFromMessage(
        serialization::SystemState const& message);
  };

  using SystemStateError = StateDifference;

  // A functor that computes f(s, y) and stores it in |variations|.
  // This functor must be called with |std::get<i>(variations).size()| equal to
  // |std::get<i>(state).size()| for all i, but there is no requirement on the
  // values in |variations|.
  RightHandSideComputation compute_derivative;
};

// A differential equation of the form X′ = A(X, t) + B(X, t), where exp(hA) and
// exp(hB) are known.  |State| is the type of X.  These equations can be solved
// using splitting methods.
template<typename... StateElements>
struct DecomposableFirstOrderDifferentialEquation final {
  static constexpr std::int64_t order = 1;
  using State = std::tuple<std::vector<StateElements>...>;

  using Flow = std::function<absl::Status(Instant const& t_initial,
                                          Instant const& t_final,
                                          State const& initial_state,
                                          State& final_state)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(Instant const& t, State const& y);

    DoublePrecision<Instant> time;
    std::tuple<std::vector<DoublePrecision<StateElements>>...> y;

    friend bool operator==(SystemState const& lhs, SystemState const& rhs) {
      return lhs.time == rhs.time && lhs.y == rhs.y;
    }
  };

  // We cannot use |Difference<StateElements>| here for the same reason.  For
  // some reason |DoublePrecision<StateElements>| above works...
  using SystemStateError =
      std::tuple<std::vector<Difference<StateElements, StateElements>>...>;

  // left_flow(t₀, t₁, X₀, X₁) sets X₁ to exp((t₁-t₀)A)X₀, and
  // right_flow(t₀, t₁, X₀, X₁) sets X₁ to exp((t₁-t₀)B)X₀.
  // The |std::vectors| in X₁ must have the same |size()| as those in X₀.  There
  // is no other requirement on their values.
  Flow left_flow;
  Flow right_flow;
};

// A differential equation of the form q″ = f(t, q, q′).
// |Position| is the type of q.
template<typename Position_>
struct ExplicitSecondOrderOrdinaryDifferentialEquation final {
  static constexpr std::int64_t order = 2;
  using IndependentVariable = Instant;
  using IndependentVariableDifference = Time;
  using Position = Position_;
  // The type of Δq.
  using Displacement = Difference<Position>;
  // The type of q′.
  using Velocity = Variation<Position>;
  // The type of q″.
  using Acceleration = Variation<Velocity>;
  using RightHandSideComputation =
      std::function<absl::Status(Instant const& t,
                                 std::vector<Position> const& positions,
                                 std::vector<Velocity> const& velocities,
                                 std::vector<Acceleration>& accelerations)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(Instant const& t,
                std::vector<Position> const& q,
                std::vector<Velocity> const& v);

    DoublePrecision<Instant> time;
    std::vector<DoublePrecision<Position>> positions;
    std::vector<DoublePrecision<Velocity>> velocities;

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

  // A functor that computes f(t, q, q′) and stores it in |accelerations|.
  // This functor must be called with |accelerations.size()| equal to
  // |positions.size()| and |velocities.size()| but there is no requirement on
  // the values in |accelerations|.
  RightHandSideComputation compute_acceleration;
};

// A differential equation of the form q″ = f(t, q).
// |Position| is the type of q.
template<typename Position_>
struct SpecialSecondOrderDifferentialEquation final {
  static constexpr std::int64_t order = 2;
  using IndependentVariable = Instant;
  using IndependentVariableDifference = Time;
  using Position = Position_;
  // The type of Δq.
  using Displacement = Difference<Position>;
  // The type of q′.
  using Velocity = Variation<Position>;
  // The type of q″.
  using Acceleration = Variation<Velocity>;
  using RightHandSideComputation =
      std::function<
          absl::Status(Instant const& t,
                       std::vector<Position> const& positions,
                       std::vector<Acceleration>& accelerations)>;

  using SystemState = typename ExplicitSecondOrderOrdinaryDifferentialEquation<
      Position>::SystemState;
  using SystemStateError =
      typename ExplicitSecondOrderOrdinaryDifferentialEquation<
          Position>::SystemStateError;

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

using internal_ordinary_differential_equations::
    DecomposableFirstOrderDifferentialEquation;
using internal_ordinary_differential_equations::
    ExplicitFirstOrderOrdinaryDifferentialEquation;
using internal_ordinary_differential_equations::
    ExplicitSecondOrderOrdinaryDifferentialEquation;
using internal_ordinary_differential_equations::IntegrationProblem;
using internal_ordinary_differential_equations::
    SpecialSecondOrderDifferentialEquation;

}  // namespace integrators
}  // namespace principia

#include "integrators/ordinary_differential_equations_body.hpp"
