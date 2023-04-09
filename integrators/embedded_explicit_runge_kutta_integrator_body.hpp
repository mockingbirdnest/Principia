#pragma once

#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "base/for_all_of.hpp"
#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace _embedded_explicit_runge_kutta_integrator {
namespace internal {

using namespace principia::base::_for_all_of;
using namespace principia::base::_not_null;
using namespace principia::geometry::_sign;
using namespace principia::numerics::_double_precision;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

template<typename Method, typename ODE_>
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::
EmbeddedExplicitRungeKuttaIntegrator() {
  // The first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance [DEP87a], equation 3.1.
    CHECK_EQ(1.0, c_[stages_ - 1]);
    CHECK_EQ(0.0, b̂_[stages_ - 1]);
    for (int j = 0; j < stages_ - 1; ++j) {
      CHECK_EQ(b̂_[j], a_(stages_ - 1, j));
    }
  }
}

template<typename Method, typename ODE_>
absl::Status
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
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
  auto const& b̂ = integrator_.b̂_;
  auto const& b = integrator_.b_;
  auto const& c = integrator_.c_;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto& first_use = this->first_use_;
  auto& parameters = this->parameters_;
  auto const& equation = this->equation_;

  // |current_state| gets updated as the integration progresses to allow
  // restartability.

  // State before the last, truncated step.
  std::optional<typename ODE::State> final_state;

  // Argument checks.
  Sign const integration_direction = Sign(parameters.first_step);
  if (integration_direction.is_positive()) {
    // Integrating forward.
    CHECK_LT(current_state.s.value, s_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.s.value, s_final);
  }
  CHECK(first_use || !parameters.last_step_is_exact)
      << "Cannot reuse an instance where the last step is exact";
  first_use = false;

  // Step.  Updated as the integration progresses to allow restartability.
  IndependentVariableDifference& h = this->step_;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<IndependentVariable>& s = current_state.s;

  // DependentVariables increment (high-order).
  DependentVariableDifferences Δŷ;
  // Current state.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  auto& ŷ = current_state.y;

  // Difference between the low- and high-order approximations.
  typename ODE::State::Error error_estimate;

  // Current Runge-Kutta stage.
  DependentVariables y_stage;

  DependentVariableDerivatives f;
  DependentVariableDerivatives last_f;
  std::vector<DependentVariableDifferences> k(stages_);

  for_all_of(ŷ, y_stage).loop([](auto const& ŷ, auto& y_stage) {
    y_stage = ŷ.value;
  });

  bool at_end = false;
  double tolerance_to_error_ratio;

  // The number of steps already performed.
  std::int64_t step_count = 0;

  absl::Status status;
  absl::Status step_status;

  if (first_same_as_last) {
    status = equation.compute_derivative(s.value, y_stage, last_f);
  }

  // No step size control on the first step.  If this instance is being
  // restarted we already have a value of |h| suitable for the next step, based
  // on the computation of |tolerance_to_error_ratio_| during the last
  // invocation.
  goto runge_kutta_step;

  while (!at_end) {
    // Compute the next step with decreasing step sizes until the error is
    // tolerable.
    do {
      // Reset the status as any error returned by a force computation for a
      // rejected step is now moot.
      step_status = absl::OkStatus();

      // Adapt step size.
      // TODO(egg): find out whether there's a smarter way to compute that root,
      // especially since we make the order compile-time.
      h *= parameters.safety_factor *
               std::pow(tolerance_to_error_ratio, 1.0 / (lower_order + 1));
      // TODO(egg): should we check whether it vanishes in double precision
      // instead?
      if (s.value + (s.error + h) == s.value) {
        return absl::Status(termination_condition::VanishingStepSize,
                            "At time " + DebugString(s.value) +
                                ", step size is effectively zero.  "
                                "Singularity or stiff system suspected.");
      }

    runge_kutta_step:
      // Termination condition.
      if (parameters.last_step_is_exact) {
        IndependentVariableDifference const s_to_end =
            (s_final - s.value) - s.error;
        at_end = integration_direction * h >=
                 integration_direction * s_to_end;
        if (at_end) {
          // The chosen step size will overshoot.  Clip it to just reach the
          // end, and terminate if the step is accepted.  Note that while this
          // step size is a good approximation, there is no guarantee that it
          // won't over/undershoot, so we still need to special case the very
          // last stage below.
          h = s_to_end;
          final_state = current_state;
        }
      }

      // Runge-Kutta iteration; fills |k|.
      for (int i = 0; i < stages_; ++i) {
        if (i == 0 && first_same_as_last) {
          // TODO(phl): Use pointers to avoid copying big objects.
          f = last_f;
        } else {
          IndependentVariable const s_stage =
              (parameters.last_step_is_exact && at_end && c[i] == 1.0)
                  ? s_final
                  : s.value + (s.error + c[i] * h);

          // TODO(phl): Should dimension |Σⱼ_aᵢⱼ_kⱼ| in the not FSAL case.
          DependentVariableDifferences Σⱼ_aᵢⱼ_kⱼ{};
          for (int j = 0; j < i; ++j) {
            for_all_of(k[j], Σⱼ_aᵢⱼ_kⱼ)
                .loop([&a, i, j](auto const& kⱼ, auto& Σⱼ_aᵢⱼ_kⱼ) {
                  Σⱼ_aᵢⱼ_kⱼ += a(i, j) * kⱼ;
                });
          }
          for_all_of(ŷ, Σⱼ_aᵢⱼ_kⱼ, y_stage)
              .loop([](auto const& ŷ, auto const& Σⱼ_aᵢⱼ_kⱼ, auto& y_stage) {
                y_stage = ŷ.value + Σⱼ_aᵢⱼ_kⱼ;
              });

          termination_condition::UpdateWithAbort(
              equation.compute_derivative(s_stage, y_stage, f), step_status);
        }
        for_all_of(f, k[i]).loop([h](auto const& f, auto& kᵢ) {
          kᵢ = h * f;
        });
      }

      // Increment computation and step size control.
      DependentVariableDifferences Σᵢ_b̂ᵢ_kᵢ{};
      DependentVariableDifferences Σᵢ_bᵢ_kᵢ{};
      for (int i = 0; i < stages_; ++i) {
        for_all_of(k[i], ŷ, Δŷ, Σᵢ_b̂ᵢ_kᵢ, Σᵢ_bᵢ_kᵢ, error_estimate)
            .loop([&b, &b̂, i](auto const& kᵢ,
                              auto const& ŷ,
                              auto& Δŷ,
                              auto& Σᵢ_b̂ᵢ_kᵢ,
                              auto& Σᵢ_bᵢ_kᵢ,
                              auto& error_estimate) {
              Σᵢ_b̂ᵢ_kᵢ += b̂[i] * kᵢ;
              Σᵢ_bᵢ_kᵢ += b[i] * kᵢ;
              Δŷ = Σᵢ_b̂ᵢ_kᵢ;
              auto const Δy = Σᵢ_bᵢ_kᵢ;
              error_estimate = Δy - Δŷ;
            });
      }
      tolerance_to_error_ratio =
          this->tolerance_to_error_ratio_(h, current_state, error_estimate);
    } while (tolerance_to_error_ratio < 1.0);

    status.Update(step_status);

    if (!parameters.last_step_is_exact && s.value + (s.error + h) > s_final) {
      // We did overshoot.  Drop the point that we just computed and exit.
      final_state = current_state;
      break;
    }

    if (first_same_as_last) {
      last_f = f;
    }

    // Increment the solution with the high-order approximation.
    s.Increment(h);
    for_all_of(Δŷ, ŷ).loop([](auto const& Δŷ, auto& ŷ) {
      ŷ.Increment(Δŷ);
    });

    RETURN_IF_STOPPED;
    append_state(current_state);
    ++step_count;
    if (absl::IsAborted(step_status)) {
      return step_status;
    } else if (step_count == parameters.max_steps && !at_end) {
      return absl::Status(termination_condition::ReachedMaximalStepCount,
                          "Reached maximum step count " +
                              std::to_string(parameters.max_steps) +
                              " at time " + DebugString(s.value) +
                              "; requested s_final is " + DebugString(s_final) +
                              ".");
    }
  }
  // The resolution is restartable from the last non-truncated state.
  CHECK(final_state);
  current_state = *final_state;
  return status;
}

template<typename Method, typename ODE_>
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_> const&
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
WriteToMessage(not_null<serialization::IntegratorInstance*> message) const {
  AdaptiveStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  [[maybe_unused]] auto* const extension =
      message
          ->MutableExtension(
              serialization::AdaptiveStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::
                  EmbeddedExplicitRungeKuttaNystromIntegratorInstance::
                      extension);
}

#if 0
template<typename Method, typename ODE_>
template<typename, typename>
not_null<std::unique_ptr<
    typename EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance>>
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
ReadFromMessage(serialization::
                    EmbeddedExplicitRungeKuttaNystromIntegratorInstance const&
                        extension,
                InitialValueProblem<ODE> const& problem,
                AppendState const& append_state,
                ToleranceToErrorRatio const& tolerance_to_error_ratio,
                Parameters const& parameters,
                Time const& time_step,
                bool const first_use,
                EmbeddedExplicitRungeKuttaIntegrator const& integrator) {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(new Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                time_step,
                                                first_use,
                                                integrator));
}
#endif

template<typename Method, typename ODE_>
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::Instance::
Instance(InitialValueProblem<ODE> const& problem,
         AppendState const& append_state,
         ToleranceToErrorRatio const& tolerance_to_error_ratio,
         Parameters const& parameters,
         typename ODE::IndependentVariableDifference const& step,
         bool const first_use,
         EmbeddedExplicitRungeKuttaIntegrator const& integrator)
    : AdaptiveStepSizeIntegrator<ODE>::Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                step,
                                                first_use),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::
NewInstance(InitialValueProblem<ODE> const& problem,
            AppendState const& append_state,
            ToleranceToErrorRatio const& tolerance_to_error_ratio,
            Parameters const& parameters) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem,
                   append_state,
                   tolerance_to_error_ratio,
                   parameters,
                   /*step=*/parameters.first_step,
                   /*first_use=*/true,
                   *this));
}

template<typename Method, typename ODE_>
void EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_>::
WriteToMessage(
    not_null<serialization::AdaptiveStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal

template<typename Method, typename ODE_>
internal::EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_> const&
EmbeddedExplicitRungeKuttaIntegrator() {
  static_assert(
      std::is_base_of<methods::EmbeddedExplicitRungeKutta, Method>::value,
      "Method must be derived from EmbeddedExplicitRungeKutta");
  static internal::EmbeddedExplicitRungeKuttaIntegrator<Method, ODE_> const
      integrator;
  return integrator;
}

}  // namespace _embedded_explicit_runge_kutta_integrator
}  // namespace integrators
}  // namespace principia
