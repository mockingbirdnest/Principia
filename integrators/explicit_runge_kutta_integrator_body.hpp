#pragma once

#include "integrators/explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_explicit_runge_kutta_integrator {

using numerics::DoublePrecision;
using namespace principia::base::_for_all_of;
using namespace principia::base::_not_null;
using namespace principia::geometry::_sign;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Method, typename ODE_>
ExplicitRungeKuttaIntegrator<Method, ODE_>::
ExplicitRungeKuttaIntegrator() {
  // The first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance [DEP87a], equation 3.1.
    CHECK_EQ(1.0, c_[stages_ - 1]);
    CHECK_EQ(0.0, b_[stages_ - 1]);
    for (int j = 0; j < stages_ - 1; ++j) {
      CHECK_EQ(b_[j], a_(stages_ - 1, j));
    }
  }
}

template<typename Method, typename ODE_>
absl::Status
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Solve(typename ODE::IndependentVariable const& s_final) {
  using IndependentVariable = typename ODE::IndependentVariable;
  using IndependentVariableDifference =
      typename ODE::IndependentVariableDifference;
  using DependentVariables = typename ODE::DependentVariables;
  using DependentVariableDifferences =
      typename ODE::DependentVariableDifferences;
  using DependentVariableDerivatives =
      typename ODE::DependentVariableDerivatives;

  auto const& a = integrator_.a_;
  auto const& b = integrator_.b_;
  auto const& c = integrator_.c_;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto const& equation = this->equation_;
  auto const& step = this->step_;

  // Argument checks.
  Sign const integration_direction = Sign(step);
  if (integration_direction.is_positive()) {
    // Integrating forward.
    CHECK_LT(current_state.s.value, s_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.s.value, s_final);
  }

  // Step.
  IndependentVariableDifference const& h = step;
  IndependentVariableDifference const abs_h = integration_direction * h;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<IndependentVariable>& s = current_state.s;

  // DependentVariables increment (high-order).
  DependentVariableDifferences Δy;
  // Current state.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  auto& y = current_state.y;

  // Current Runge-Kutta stage.
  DependentVariables y_stage;

  DependentVariableDerivatives f;
  DependentVariableDerivatives last_f;
  std::vector<DependentVariableDifferences> k(stages_);

  for_all_of(y, y_stage).loop([](auto const& y, auto& y_stage) {
    y_stage = y.value;
  });

  absl::Status status;

  if (first_same_as_last) {
    status = equation.compute_derivative(s.value, y_stage, last_f);
  }

  while (abs_h <= Abs((s_final - s.value) - s.error)) {
    // Runge-Kutta iteration; fills |k|.
    for (int i = 0; i < stages_; ++i) {
      if (i == 0 && first_same_as_last) {
        // TODO(phl): Use pointers to avoid copying big objects.
        f = last_f;
      } else {
        DependentVariableDifferences Σⱼ_aᵢⱼ_kⱼ{};

        for (int j = 0; j < i; ++j) {
          for_all_of(k[j], y, y_stage, Σⱼ_aᵢⱼ_kⱼ)
              .loop([&a, i, j](auto const& kⱼ,
                               auto const& y,
                               auto& y_stage,
                               auto& Σⱼ_aᵢⱼ_kⱼ) {
                Σⱼ_aᵢⱼ_kⱼ += a(i, j) * kⱼ;
              });
        }
        for_all_of(y, Σⱼ_aᵢⱼ_kⱼ, y_stage)
            .loop([](auto const& y, auto const& Σⱼ_aᵢⱼ_kⱼ, auto& y_stage) {
              y_stage = y.value + Σⱼ_aᵢⱼ_kⱼ;
            });

        termination_condition::UpdateWithAbort(
            equation.compute_derivative(
                s.value + (s.error + c[i] * h), y_stage, f),
            status);
      }
      for_all_of(f, k[i]).loop([h](auto const& f, auto& kᵢ) {
        kᵢ = h * f;
      });
    }

    // Increment computation.
    DependentVariableDifferences Σᵢ_bᵢ_kᵢ{};
    for (int i = 0; i < stages_; ++i) {
      for_all_of(k[i], y, Δy, Σᵢ_bᵢ_kᵢ)
          .loop([&a, &b, i](
                    auto const& kᵢ, auto const& y, auto& Δy, auto& Σᵢ_bᵢ_kᵢ) {
            Σᵢ_bᵢ_kᵢ += b[i] * kᵢ;
            Δy = Σᵢ_bᵢ_kᵢ;
          });
    }

    if (first_same_as_last) {
      last_f = f;
    }

    // Increment the solution.
    s.Increment(h);
    for_all_of(Δy, y).loop([](auto const& Δy, auto& y) {
      y.Increment(Δy);
    });

    RETURN_IF_STOPPED;
    append_state(current_state);
    if (absl::IsAborted(status)) {
      return status;
    }
  }

  return status;
}

template<typename Method, typename ODE_>
ExplicitRungeKuttaIntegrator<Method, ODE_> const&
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
ExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Instance(InitialValueProblem<ODE> const& problem,
         AppendState const& append_state,
         typename ODE::IndependentVariableDifference const& step,
         ExplicitRungeKuttaIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem,
                                             append_state,
                                             step),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
ExplicitRungeKuttaIntegrator<Method, ODE_>::
NewInstance(InitialValueProblem<ODE> const& problem,
            AppendState const& append_state,
            typename ODE::IndependentVariableDifference const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename ODE_>
void ExplicitRungeKuttaIntegrator<Method, ODE_>::
WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal_explicit_runge_kutta_integrator

template<typename Method, typename ODE_>
internal_explicit_runge_kutta_integrator::
    ExplicitRungeKuttaIntegrator<Method, ODE_> const&
ExplicitRungeKuttaIntegrator() {
  static_assert(
      std::is_base_of<methods::ExplicitRungeKutta,
                      Method>::value,
      "Method must be derived from ExplicitRungeKutta");
  static internal_explicit_runge_kutta_integrator::
      ExplicitRungeKuttaIntegrator<Method, ODE_> const integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
