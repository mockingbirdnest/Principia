#pragma once

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <vector>

#include "base/jthread.hpp"
#include "geometry/sign.hpp"
#include "integrators/methods.hpp"
#include "numerics/ulp_distance.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using numerics::DoublePrecision;
using numerics::ULPDistance;
using quantities::Abs;
using namespace principia::base::_not_null;
using namespace principia::geometry::_sign;

template<typename Method, typename ODE_>
absl::Status SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::Solve(Instant const& t_final) {
  using Position = typename ODE::DependentVariable;
  using Displacement = typename ODE::DependentVariableDifference;
  using Velocity = typename ODE::DependentVariableDerivative;
  using Acceleration = typename ODE::DependentVariableDerivative2;

  auto const& a = integrator_.a_;
  auto const& b = integrator_.b_;
  auto const& c = integrator_.c_;

  auto& current_state = this->current_state_;
  auto& append_state = this->append_state_;
  auto const& equation = this->equation_;
  auto const& step = this->step_;

  // |current_state| is updated as the integration progresses to allow
  // restartability.

  // Argument checks.
  int const dimension = current_state.positions.size();
  CHECK_NE(Time(), step);
  Sign const integration_direction = Sign(step);
  if (integration_direction.is_positive()) {
    // Integrating forward.
    CHECK_LT(current_state.time.value, t_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.time.value, t_final);
  }

  // Time step.
  Time const& h = step;
  Time const abs_h = integration_direction * h;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment.
  std::vector<Displacement> Δq(dimension);
  // Velocity increment.
  std::vector<Velocity> Δv(dimension);
  // Current position.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Position>>& q = current_state.positions;
  // Current velocity.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Velocity>>& v = current_state.velocities;

  // Current Runge-Kutta-Nyström stage.
  std::vector<Position> q_stage(dimension);
  // Accelerations at the current stage.
  std::vector<Acceleration> g(dimension);

  // The first full stage of the step, i.e. the first stage where
  // exp(bᵢ h B) exp(aᵢ h A) must be entirely computed.
  // 0 in the non-FSAL BA case, 1 in the ABA case since b₀ = 0 means the first
  // stage is only exp(a₀ h A), 1 in the BAB case, since the previous
  // right-hand-side evaluation can be used for exp(bᵢ h B).  Note that in the
  // BAB case, we need to start things with an evaluation since there is no
  // previous evaluation.
  constexpr int first_stage = composition == BA ? 0 : 1;

  absl::Status status;

  if (composition == BAB) {
    for (int k = 0; k < dimension; ++k) {
      q_stage[k] = q[k].value;
    }
    status.Update(equation.compute_acceleration(t.value, q_stage, g));
  }

  while (abs_h <= Abs((t_final - t.value) - t.error)) {
    std::fill(Δq.begin(), Δq.end(), Displacement{});
    std::fill(Δv.begin(), Δv.end(), Velocity{});

    if (first_stage == 1) {
      for (int k = 0; k < dimension; ++k) {
        if (composition == BAB) {
          // exp(b₀ h B)
          Δv[k] += h * b[0] * g[k];
        }
        // exp(a₀ h A)
        Δq[k] += h * a[0] * (v[k].value + Δv[k]);
      }
    }

    for (int i = first_stage; i < stages_; ++i) {
      for (int k = 0; k < dimension; ++k) {
        q_stage[k] = q[k].value + Δq[k];
      }
      termination_condition::UpdateWithAbort(
          equation.compute_acceleration(
              t.value + (t.error + c[i] * h), q_stage, g),
          status);
      for (int k = 0; k < dimension; ++k) {
        // exp(bᵢ h B)
        Δv[k] += h * b[i] * g[k];
        // NOTE(egg): in the BAB case, at the last stage, this will be an
        // exercise in adding 0.  I don't think the optimizer can know that.  Do
        // we care?
        // exp(aᵢ h A)
        Δq[k] += h * a[i] * (v[k].value + Δv[k]);
      }
    }

    // Increment the solution.
    t.Increment(h);
    for (int k = 0; k < dimension; ++k) {
      q[k].Increment(Δq[k]);
      v[k].Increment(Δv[k]);
    }
    RETURN_IF_STOPPED;
    append_state(current_state);
    if (absl::IsAborted(status)) {
      return status;
    }
  }

  return status;
}

template<typename Method, typename ODE_>
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_> const&
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename ODE_>
void SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  [[maybe_unused]] auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::SymplecticRungeKuttaNystromIntegratorInstance::
                  extension);
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<
    typename SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::Instance>>
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::Instance::
ReadFromMessage(
    serialization::SymplecticRungeKuttaNystromIntegratorInstance const&
        extension,
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    SymplecticRungeKuttaNyströmIntegrator const& integrator) {
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, integrator));
}

template<typename Method, typename ODE_>
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::
Instance::Instance(InitialValueProblem<ODE> const& problem,
                   AppendState const& append_state,
                   Time const& step,
                   SymplecticRungeKuttaNyströmIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem,
                                             std::move(append_state),
                                             step),
      integrator_(integrator) {}

template<typename Method, typename ODE_>
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::
SymplecticRungeKuttaNyströmIntegrator() {
  DoublePrecision<double> c_i(0.0);
  for (int i = 0; i < stages_; ++i) {
    c_[i] = c_i.value;
    c_i += DoublePrecision<double>(a_[i]);
  }
  CHECK_LE(ULPDistance(1.0, c_i.value), 4);
  if (composition == ABA) {
    CHECK_EQ(0.0, b_[0]);
  } else if (composition == BAB) {
    CHECK_EQ(0.0, a_[stages_ - 1]);
  }
  if (time_reversible) {
    switch (composition) {
      case ABA:
        for (int i = 0; i < stages_; ++i) {
          CHECK_EQ(a_[i], a_[stages_ - 1 - i]);
        }
        for (int i = 0; i < stages_ - 1; ++i) {
          CHECK_EQ(b_[i + 1], b_[stages_ - 1 - i]);
        }
        break;
      case BAB:
        for (int i = 0; i < stages_ - 1; ++i) {
          CHECK_EQ(a_[i], a_[stages_ - 2 - i]);
        }
        for (int i = 0; i < stages_; ++i) {
          CHECK_EQ(b_[i], b_[stages_ - 1 - i]);
        }
        break;
      case BA:
        LOG(FATAL) << "Time-reversible compositions have the FSAL property";
        break;
      default:
        LOG(FATAL) << "Invalid CompositionMethod";
    }
  }
}

template<typename Method, typename ODE_>
not_null<std::unique_ptr<typename Integrator<ODE_>::Instance>>
SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::NewInstance(
    InitialValueProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename ODE_>
void SymplecticRungeKuttaNyströmIntegrator<Method, ODE_>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
  message->set_composition_method(composition);
}

}  // namespace internal_symplectic_runge_kutta_nyström_integrator

template<typename Method, typename ODE_>
internal_symplectic_runge_kutta_nyström_integrator::
    SymplecticRungeKuttaNyströmIntegrator<Method, ODE_> const&
SymplecticRungeKuttaNyströmIntegrator() {
  static_assert(
      std::is_base_of<methods::SymplecticRungeKuttaNyström, Method>::value,
      "Method must be derived from SymplecticRungeKuttaNyström");
  static internal_symplectic_runge_kutta_nyström_integrator::
      SymplecticRungeKuttaNyströmIntegrator<Method, ODE_> const integrator;
  return integrator;
}

template<typename Method,
         serialization::FixedStepSizeIntegrator::CompositionMethod composition,
         typename ODE_>
internal_symplectic_runge_kutta_nyström_integrator::
    SymplecticRungeKuttaNyströmIntegrator<
        typename methods::AsSymplecticRungeKuttaNyström<Method,
                                                        composition>::Method,
        ODE_> const&
SymplecticRungeKuttaNyströmIntegrator() {
  static_assert(
      std::is_base_of<methods::SymplecticPartitionedRungeKutta, Method>::value,
      "Method must be derived from SymplecticPartitionedRungeKutta");
  static internal_symplectic_runge_kutta_nyström_integrator::
      SymplecticRungeKuttaNyströmIntegrator<
          typename methods::AsSymplecticRungeKuttaNyström<Method,
                                                          composition>::Method,
          ODE_> const integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
