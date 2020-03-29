
#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include <algorithm>
#include <list>
#include <utility>
#include <vector>

#include "geometry/serialization.hpp"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {
namespace internal_symmetric_linear_multistep_integrator {

using base::make_not_null_unique;
using geometry::QuantityOrMultivectorSerializer;

int const startup_step_divisor = 16;

template<typename Method, typename Position>
Status SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Solve(
    Instant const& t_final) {
  using Acceleration = typename ODE::Acceleration;
  using Displacement = typename ODE::Displacement;
  using DoubleDisplacement = DoublePrecision<Displacement>;
  using DoubleDisplacements = std::vector<DoubleDisplacement>;
  using DoublePosition = DoublePrecision<Position>;

  auto const& ɑ = integrator_.ɑ_;
  auto const& β_numerator = integrator_.β_numerator_;
  auto const& β_denominator = integrator_.β_denominator_;

  auto& current_state = this->current_state_;
  auto& append_state = this->append_state_;
  auto const& step = this->step_;
  auto const& equation = this->equation_;

  if (previous_steps_.size() < order) {
    StartupSolve(t_final);

    // If |t_final| is not large enough, we may not have generated enough
    // points.  Bail out, we'll continue the next time |Solve| is called.
    if (previous_steps_.size() < order) {
      return Status::OK;
    }
  }
  CHECK_EQ(previous_steps_.size(), order);

  // Argument checks.
  int const dimension = previous_steps_.back().displacements.size();

  // Time step.
  CHECK_LT(Time(), step);
  Time const& h = step;
  // Current time.
  DoublePrecision<Instant> t = previous_steps_.back().time;
  // Order.
  int const k = order;

  Status status;
  std::vector<Position> positions(dimension);

  DoubleDisplacements Σj_minus_ɑj_qj(dimension);
  std::vector<Acceleration> Σj_βj_numerator_aj(dimension);
  while (h <= (t_final - t.value) - t.error) {
    // We take advantage of the symmetry to iterate on the list of previous
    // steps from both ends.
    auto front_it = previous_steps_.begin();
    auto back_it = previous_steps_.rbegin();

    // This block corresponds to j = 0.  We must not pair it with j = k.
    {
      DoubleDisplacements const& qj = front_it->displacements;
      std::vector<Acceleration> const& aj = front_it->accelerations;
      double const ɑj = ɑ[0];
      double const βj_numerator = β_numerator[0];
      for (int d = 0; d < dimension; ++d) {
        Σj_minus_ɑj_qj[d] = Scale(-ɑj, qj[d]);
        Σj_βj_numerator_aj[d] = βj_numerator * aj[d];
      }
      ++front_it;
    }
    // The generic value of j, paired with k - j.
    for (int j = 1; j < k / 2; ++j) {
      DoubleDisplacements const& qj = front_it->displacements;
      DoubleDisplacements const& qk_minus_j = back_it->displacements;
      std::vector<Acceleration> const& aj = front_it->accelerations;
      std::vector<Acceleration> const& ak_minus_j = back_it->accelerations;
      double const ɑj = ɑ[j];
      double const βj_numerator = β_numerator[j];
      for (int d = 0; d < dimension; ++d) {
        Σj_minus_ɑj_qj[d] -= Scale(ɑj, qj[d]);
        Σj_minus_ɑj_qj[d] -= Scale(ɑj, qk_minus_j[d]);
        Σj_βj_numerator_aj[d] += βj_numerator * (aj[d] + ak_minus_j[d]);
      }
      ++front_it;
      ++back_it;
    }
    // This block corresponds to j = k / 2.  We must not pair it with j = k / 2.
    {
      DoubleDisplacements const& qj = front_it->displacements;
      std::vector<Acceleration> const& aj = front_it->accelerations;
      double const ɑj = ɑ[k / 2];
      double const βj_numerator = β_numerator[k / 2];
      for (int d = 0; d < dimension; ++d) {
        Σj_minus_ɑj_qj[d] -= Scale(ɑj, qj[d]);
        Σj_βj_numerator_aj[d] += βj_numerator * aj[d];
      }
    }

    // Create a new step in the instance.
    t.Increment(h);
    previous_steps_.emplace_back();
    Step& current_step = previous_steps_.back();
    current_step.time = t;
    current_step.accelerations.resize(dimension);
    current_step.displacements.reserve(dimension);

    // Fill the new step.  We skip the division by ɑk as it is equal to 1.0.
    double const ɑk = ɑ[0];
    DCHECK_EQ(ɑk, 1.0);
    for (int d = 0; d < dimension; ++d) {
      DoubleDisplacement& current_displacement = Σj_minus_ɑj_qj[d];
      current_displacement.Increment(h * h *
                                     Σj_βj_numerator_aj[d] / β_denominator);
      current_step.displacements.push_back(current_displacement);
      DoublePosition const current_position =
          DoublePosition() + current_displacement;
      positions[d] = current_position.value;
      current_state.positions[d] = current_position;
    }
    status.Update(equation.compute_acceleration(t.value,
                                                positions,
                                                current_step.accelerations));
    previous_steps_.pop_front();

    ComputeVelocityUsingCohenHubbardOesterwinter();

    // Inform the caller of the new state.
    current_state.time = t;
    append_state(current_state);
  }

  return status;
}

template<typename Method, typename Position>
SymmetricLinearMultistepIntegrator<Method, Position> const&
SymmetricLinearMultistepIntegrator<Method, Position>::Instance::integrator()
    const {
  return integrator_;
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename Position>
void SymmetricLinearMultistepIntegrator<Method, Position>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::SymmetricLinearMultistepIntegratorInstance::
                  extension);
  for (auto const& previous_step : previous_steps_) {
    previous_step.WriteToMessage(extension->add_previous_steps());
  }
  extension->set_startup_step_index(startup_step_index_);
}

template<typename Method, typename Position>
template<typename, typename>
not_null<std::unique_ptr<
    typename SymmetricLinearMultistepIntegrator<Method, Position>::Instance>>
SymmetricLinearMultistepIntegrator<Method, Position>::Instance::ReadFromMessage(
    serialization::SymmetricLinearMultistepIntegratorInstance const& extension,
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    SymmetricLinearMultistepIntegrator const& integrator) {
  std::list<Step> previous_steps;
  for (auto const& previous_step : extension.previous_steps()) {
    previous_steps.push_back(Step::ReadFromMessage(previous_step));
  }
  return std::unique_ptr<Instance>(new Instance(problem,
                                                append_state,
                                                step,
                                                extension.startup_step_index(),
                                                previous_steps,
                                                integrator));
}

template<typename Method, typename Position>
void SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Step::
WriteToMessage(
    not_null<serialization::SymmetricLinearMultistepIntegratorInstance::
                  Step*> const message) const {
  using AccelerationSerializer = QuantityOrMultivectorSerializer<
      typename ODE::Acceleration,
      serialization::SymmetricLinearMultistepIntegratorInstance::Step::
          Acceleration>;
  for (auto const& displacement : displacements) {
    displacement.WriteToMessage(message->add_displacements());
  }
  for (auto const& acceleration : accelerations) {
    AccelerationSerializer::WriteToMessage(acceleration,
                                           message->add_accelerations());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename Method, typename Position>
template<typename, typename>
typename SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Step
SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Step::
ReadFromMessage(
    serialization::SymmetricLinearMultistepIntegratorInstance::Step const&
        message) {
  using AccelerationSerializer = QuantityOrMultivectorSerializer<
      typename ODE::Acceleration,
      serialization::SymmetricLinearMultistepIntegratorInstance::Step::
          Acceleration>;
  Step step;
  for (auto const& displacement : message.displacements()) {
    step.displacements.push_back(
        DoublePrecision<typename ODE::Displacement>::ReadFromMessage(
            displacement));
  }
  for (auto const& acceleration : message.accelerations()) {
    step.accelerations.push_back(
        AccelerationSerializer::ReadFromMessage(acceleration));
  }
  step.time = DoublePrecision<Instant>::ReadFromMessage(message.time());
  return step;
}

template<typename Method, typename Position>
SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    SymmetricLinearMultistepIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      integrator_(integrator) {
  previous_steps_.emplace_back();
  FillStepFromSystemState(this->equation_,
                          this->current_state_,
                          this->previous_steps_.back());
}

template<typename Method, typename Position>
SymmetricLinearMultistepIntegrator<Method, Position>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    int const startup_step_index,
    std::list<Step> previous_steps,
    SymmetricLinearMultistepIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem, append_state, step),
      startup_step_index_(startup_step_index),
      previous_steps_(std::move(previous_steps)),
      integrator_(integrator) {}

template<typename Method, typename Position>
void SymmetricLinearMultistepIntegrator<Method, Position>::
Instance::StartupSolve(Instant const& t_final) {
  auto& current_state = this->current_state_;
  auto const& step = this->step_;
  auto const& equation = this->equation_;

  Time const startup_step = step / startup_step_divisor;

  CHECK(!previous_steps_.empty());
  CHECK_LT(previous_steps_.size(), order);

  auto const startup_append_state =
      [this](typename ODE::SystemState const& state) {
        // Stop changing anything once we're done with the startup.  We may be
        // called one more time by the |startup_integrator_|.
        if (previous_steps_.size() < order) {
          this->current_state_ = state;
          // The startup integrator has a smaller step.  We do not record all
          // the states it computes, but only those that are a multiple of the
          // main integrator step.
          if (++startup_step_index_ % startup_step_divisor == 0) {
            CHECK_LT(previous_steps_.size(), order);
            previous_steps_.emplace_back();
            FillStepFromSystemState(this->equation_,
                                    this->current_state_,
                                    previous_steps_.back());
            // This call must happen last for a subtle reason: the callback may
            // want to |Clone| this instance (see |Ephemeris::Checkpoint|) in
            // which cases it is necessary that all the member variables be
            // filled for restartability to work.
            this->append_state_(state);
          }
        }
      };

  auto const startup_instance =
      integrator_.startup_integrator_.NewInstance({equation, current_state},
                                                  startup_append_state,
                                                  startup_step);

  startup_instance->Solve(
      std::min(current_state.time.value +
                   (order - previous_steps_.size()) * step + step / 2.0,
               t_final));

  CHECK_LE(previous_steps_.size(), order);
}

template<typename Method, typename Position>
void SymmetricLinearMultistepIntegrator<Method, Position>::
Instance::ComputeVelocityUsingCohenHubbardOesterwinter() {
  using Acceleration = typename ODE::Acceleration;
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;

  // For the computation of the velocity we use a formula similar to that of
  // Cohen, Hubbard, Oesterwinter (1973), Astronomical papers prepared for the
  // use of the American ephemeris and nautical almanac, Volume XXII, Part I.
  // More specifically, we use the coefficients η from
  // cohen_hubbard_oesterwinter.wl.
  auto const& cohen_hubbard_oesterwinter =
      integrator_.cohen_hubbard_oesterwinter_;

  int const dimension = previous_steps_.back().displacements.size();
  auto& current_state = this->current_state_;
  auto const& step = this->step_;

  current_state.velocities.reserve(dimension);
  for (int d = 0; d < dimension; ++d) {
    DoublePrecision<Velocity>& velocity = current_state.velocities[d];
    auto it = previous_steps_.rbegin();

    // Compute the displacement difference using double precision.
    DoublePrecision<Displacement> displacement_change =
        it->displacements[d] - std::next(it)->displacements[d];
    velocity = DoublePrecision<Velocity>(
        (displacement_change.value + displacement_change.error) / step);

    Acceleration weighted_accelerations;
    for (int i = 0; i < cohen_hubbard_oesterwinter.numerators.size; ++i, ++it) {
      weighted_accelerations +=
          cohen_hubbard_oesterwinter.numerators[i] * it->accelerations[d];
    }

    velocity.value +=
        weighted_accelerations * step / cohen_hubbard_oesterwinter.denominator;
  }
}

template<typename Method, typename Position>
void SymmetricLinearMultistepIntegrator<Method, Position>::
Instance::FillStepFromSystemState(ODE const& equation,
                                  typename ODE::SystemState const& state,
                                  Step& step) {
  std::vector<typename ODE::Position> positions;
  step.time = state.time;
  for (auto const& position : state.positions) {
    step.displacements.push_back(position - DoublePrecision<Position>());
    positions.push_back(position.value);
  }
  step.accelerations.resize(step.displacements.size());
  // Ignore the status here.  We are merely computing the acceleration to store
  // it, not to advance an integrator.
  equation.compute_acceleration(step.time.value,
                                positions,
                                step.accelerations);
}

template<typename Method, typename Position>
SymmetricLinearMultistepIntegrator<Method, Position>::
SymmetricLinearMultistepIntegrator(
    FixedStepSizeIntegrator<ODE> const& startup_integrator)
    : startup_integrator_(startup_integrator),
      cohen_hubbard_oesterwinter_(CohenHubbardOesterwinterOrder<order>()) {
  CHECK_EQ(ɑ_[0], 1.0);
  CHECK_EQ(β_numerator_[0], 0.0);
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
SymmetricLinearMultistepIntegrator<Method, Position>::NewInstance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Method, typename Position>
void SymmetricLinearMultistepIntegrator<Method, Position>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> message) const {
  message->set_kind(Method::kind);
}

}  // namespace internal_symmetric_linear_multistep_integrator

template<typename Method, typename Position>
internal_symmetric_linear_multistep_integrator::
    SymmetricLinearMultistepIntegrator<Method, Position> const&
    SymmetricLinearMultistepIntegrator() {
  static_assert(
      std::is_base_of<methods::SymmetricLinearMultistep, Method>::value,
      "Method must be derived from SymmetricLinearMultistep");
  // TODO(phl): Someday, and that day may never come, I will call upon you to
  // expose the startup integrator to the clients.  But until that day, accept
  // this Blanes-Moan integrator as a gift.
  static internal_symmetric_linear_multistep_integrator::
      SymmetricLinearMultistepIntegrator<Method, Position> const integrator(
          SymplecticRungeKuttaNyströmIntegrator<methods::BlanesMoan2002SRKN14A,
                                                Position>());
  return integrator;
}

}  // namespace integrators
}  // namespace principia

#undef PRINCIPIA_USE_COHEN_HUBBARD_OESTERWINTER
