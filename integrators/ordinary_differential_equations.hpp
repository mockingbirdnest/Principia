
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
using quantities::Quotient;
using quantities::Time;
using quantities::Variation;

// A differential equation of the form y′ = f(y, t).
// |State| is the type of y.
template<typename... StateTypes>
struct ExplicitFirstOrderOrdinaryDifferentialEquation final {
  using State = std::tuple<std::vector<StateTypes>...>;
  // NOTE(eggrobin): we cannot use |Variation| here because of a bug in MSVC
  // when parameter packs interact with defaulted template parameters (and
  // possibly decltype), but this is just
  //   std::tuple<std::vector<Variation<StateTypes>>...>.
  using StateVariation = std::tuple<
      std::vector<Quotient<Difference<StateTypes, StateTypes>, Time>>...>;

  using RightHandSideComputation =
      std::function<Status(Instant const& t,
                           State const& state,
                           StateVariation& derivatives)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(State const& y, Instant const& t);

    std::tuple<std::vector<DoublePrecision<StateTypes>>...> y;
    DoublePrecision<Instant> time;

    friend bool operator==(SystemState const& lhs, SystemState const& rhs) {
      return lhs.y == rhs.y && lhs.time == rhs.time;
    }
  };

  // We cannot use |Difference<StateTypes>| here for the same reason.  For some
  // reason |DoublePrecision<StateTypes>| above works...
  using SystemStateError =
      std::tuple<std::vector<Difference<StateTypes, StateTypes>>...>;

  // A functor that computes f(y, t) and stores it in |derivatives|.
  // This functor must be called with |std::get<i>(derivatives).size()| equal to
  // |std::get<i>(state).size()| for all i, but there is no requirement on the
  // values in |derivatives|.
  RightHandSideComputation compute_derivative;
};

// A differential equation of the form X′ = A(X, t) + B(X, t), where exp(hA) and
// exp(hB) are known.  |State| is the type of X.
template<typename... StateTypes>
struct Splitting final {
  using State = std::tuple<std::vector<StateTypes>...>;

  using Flow = std::function<Status(Instant const& t_initial,
                                    Instant const& t_final,
                                    State const& initial_state,
                                    State& final_state)>;

  struct SystemState final {
    SystemState() = default;
    SystemState(State const& y, Instant const& t);

    std::tuple<std::vector<DoublePrecision<StateTypes>>...> y;
    DoublePrecision<Instant> time;

    friend bool operator==(SystemState const& lhs, SystemState const& rhs) {
      return lhs.y == rhs.y && lhs.time == rhs.time;
    }
  };

  // We cannot use |Difference<StateTypes>| here for the same reason.  For some
  // reason |DoublePrecision<StateTypes>| above works...
  using SystemStateError =
      std::tuple<std::vector<Difference<StateTypes, StateTypes>>...>;

  // left_flow(t₀, t₁, X₀, X₁) sets X₁ to exp((t₁-t₀)A)X₀, and
  // right_flow(t₀, t₁, X₀, X₁) sets X₁ to exp((t₁-t₀)B)X₀.
  // The |std::vectors| in X₁ must have the same |size()| as those in X₀.  There
  // is no other requirement on their values.
  Flow left_flow;
  Flow right_flow;
};

// A differential equation of the form q″ = f(q, q′, t).
// |Position| is the type of q.
template<typename Position_>
struct ExplicitSecondOrderOrdinaryDifferentialEquation final {
  using Position = Position_;
  // The type of Δq.
  using Displacement = Difference<Position>;
  // The type of q′.
  using Velocity = Variation<Position>;
  // The type of q″.
  using Acceleration = Variation<Velocity>;
  using RightHandSideComputation =
      std::function<Status(Instant const& t,
                           std::vector<Position> const& positions,
                           std::vector<Velocity> const& velocities,
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
