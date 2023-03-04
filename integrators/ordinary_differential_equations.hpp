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

using namespace principia::base::_not_null;
using namespace principia::geometry::_named_quantities;
using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

// A differential equation of the form y′ = f(s, y).
// |DependentVariable| are the types of the elements of y.
template<typename IndependentVariable_, typename... DependentVariable>
struct ExplicitFirstOrderOrdinaryDifferentialEquation final {
  static constexpr std::int64_t order = 1;
  using IndependentVariable = IndependentVariable_;
  using IndependentVariableDifference = Difference<IndependentVariable>;
  using DependentVariables = std::tuple<DependentVariable...>;
  using DependentVariableDifferences =
      std::tuple<Difference<DependentVariable>...>;
  using DependentVariableDerivatives = std::tuple<
      Derivative<DependentVariable, IndependentVariable>...>;

  // A functor that computes f(s, y) and stores it in |y′|.  This functor must
  // be called with |std::get<i>(y′).size()| equal to |std::get<i>(y).size()|
  // for all i, but there is no requirement on the values in |y′|.
  using RightHandSideComputation =
      std::function<absl::Status(IndependentVariable const& s,
                                 DependentVariables const& y,
                                 DependentVariableDerivatives& yʹ)>;

  struct State final {
    using Error = DependentVariableDifferences;

    State() = default;
    State(IndependentVariable const& s, DependentVariables const& y);

    DoublePrecision<IndependentVariable> s;
    std::tuple<DoublePrecision<DependentVariable>...> y;

    friend bool operator==(State const& lhs, State const& rhs) {
      return lhs.y == rhs.y && lhs.s == rhs.s;
    }

    void WriteToMessage(not_null<serialization::State*> message) const;
    static State ReadFromMessage(serialization::State const& message);
  };

  RightHandSideComputation compute_derivative;
};

// A differential equation of the form X′ = A(X, t) + B(X, t), where exp(hA) and
// exp(hB) are known.  |DependentVariable| are the types of the elements of X.
// These equations can be solved using splitting methods.
template<typename... DependentVariable>
struct DecomposableFirstOrderDifferentialEquation final {
  static constexpr std::int64_t order = 1;
  using IndependentVariable = Instant;
  using IndependentVariableDifference = Time;
  using DependentVariables = std::tuple<std::vector<DependentVariable>...>;
  using DependentVariableDifferences =
      std::tuple<std::vector<Difference<DependentVariable>>...>;

  using Flow =
      std::function<absl::Status(IndependentVariable const& t_initial,
                                 IndependentVariable const& t_final,
                                 DependentVariables const& initial_state,
                                 DependentVariables& final_state)>;

  struct State final {
    using Error = DependentVariableDifferences;

    State() = default;
    State(IndependentVariable const& t, DependentVariables const& y);

    DoublePrecision<Instant> time;
    std::tuple<std::vector<DoublePrecision<DependentVariable>>...> y;

    friend bool operator==(State const& lhs, State const& rhs) {
      return lhs.time == rhs.time && lhs.y == rhs.y;
    }
  };

  // left_flow(t₀, t₁, X₀, X₁) sets X₁ to exp((t₁-t₀)A)X₀, and
  // right_flow(t₀, t₁, X₀, X₁) sets X₁ to exp((t₁-t₀)B)X₀.
  // The |std::vectors| in X₁ must have the same |size()| as those in X₀.  There
  // is no other requirement on their values.
  Flow left_flow;
  Flow right_flow;
};

// A differential equation of the form q″ = f(t, q, q′).
// |DependentVariable_| is the type of q.
template<typename DependentVariable_>
struct ExplicitSecondOrderOrdinaryDifferentialEquation final {
  static constexpr std::int64_t order = 2;
  using IndependentVariable = Instant;
  using IndependentVariableDifference = Time;
  using DependentVariable = DependentVariable_;
  // The type of Δq.
  using DependentVariableDifference = Difference<DependentVariable>;
  // The type of q′.
  using DependentVariableDerivative =
      Derivative<DependentVariable, IndependentVariable>;
  // The type of q″.
  using DependentVariableDerivative2 =
      Derivative<DependentVariable, IndependentVariable, 2>;
  using DependentVariables = std::vector<DependentVariable>;
  using DependentVariableDifferences = std::vector<DependentVariableDifference>;
  using DependentVariableDerivatives = std::vector<DependentVariableDerivative>;
  using DependentVariableDerivatives2 =
      std::vector<DependentVariableDerivative2>;

  // A functor that computes f(t, q, q′) and stores it in |accelerations|.
  // This functor must be called with |accelerations.size()| equal to
  // |positions.size()| and |velocities.size()| but there is no requirement on
  // the values in |accelerations|.
  using RightHandSideComputation =
      std::function<absl::Status(IndependentVariable const& t,
                                 DependentVariables const& positions,
                                 DependentVariableDerivatives const& velocities,
                                 DependentVariableDerivatives2& accelerations)>;

  struct State final {
    struct Error final {
      DependentVariableDifferences position_error;
      DependentVariableDerivatives velocity_error;
    };

    State() = default;
    State(IndependentVariable const& t,
          DependentVariables const& q,
          DependentVariableDerivatives const& v);

    DoublePrecision<IndependentVariable> time;
    std::vector<DoublePrecision<DependentVariable>> positions;
    std::vector<
        DoublePrecision<Derivative<DependentVariable, IndependentVariable>>>
        velocities;

    friend bool operator==(State const& lhs, State const& rhs) {
      return lhs.positions == rhs.positions &&
             lhs.velocities == rhs.velocities &&
             lhs.time == rhs.time;
    }

    void WriteToMessage(not_null<serialization::State*> message) const;
    static State ReadFromMessage(serialization::State const& message);
  };

  RightHandSideComputation compute_acceleration;
};

// A differential equation of the form q″ = f(t, q).
// |DependentVariable_| is the type of q.
template<typename DependentVariable_>
struct SpecialSecondOrderDifferentialEquation final {
  static constexpr std::int64_t order = 2;
  using IndependentVariable = Instant;
  using IndependentVariableDifference = Time;
  using DependentVariable = DependentVariable_;
  // The type of Δq.
  using DependentVariableDifference = Difference<DependentVariable>;
  // The type of q′.
  using DependentVariableDerivative =
      Derivative<DependentVariable, IndependentVariable>;
  // The type of qʺ.
  using DependentVariableDerivative2 =
      Derivative<DependentVariable, IndependentVariable, 2>;
  using DependentVariables = std::vector<DependentVariable>;
  using DependentVariableDifferences = std::vector<DependentVariableDifference>;
  using DependentVariableDerivatives2 =
      std::vector<DependentVariableDerivative2>;

  using RightHandSideComputation =
      std::function<absl::Status(IndependentVariable const& t,
                                 DependentVariables const& positions,
                                 DependentVariableDerivatives2& accelerations)>;

  using State = typename ExplicitSecondOrderOrdinaryDifferentialEquation<
      DependentVariable>::State;

  // A functor that computes f(q, t) and stores it in |accelerations|.
  // This functor must be called with |accelerations.size()| equal to
  // |positions.size()|, but there is no requirement on the values in
  // |acceleration|.
  RightHandSideComputation compute_acceleration;
};

// An initial value problem.
template<typename ODE>
struct InitialValueProblem final {
  ODE equation;
  typename ODE::State initial_state;
};

}  // namespace internal_ordinary_differential_equations

using internal_ordinary_differential_equations::
    DecomposableFirstOrderDifferentialEquation;
using internal_ordinary_differential_equations::
    ExplicitFirstOrderOrdinaryDifferentialEquation;
using internal_ordinary_differential_equations::
    ExplicitSecondOrderOrdinaryDifferentialEquation;
using internal_ordinary_differential_equations::InitialValueProblem;
using internal_ordinary_differential_equations::
    SpecialSecondOrderDifferentialEquation;

}  // namespace integrators
}  // namespace principia

#include "integrators/ordinary_differential_equations_body.hpp"
