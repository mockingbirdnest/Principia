#pragma once

#include "integrators/symmetric_linear_multistep_integrator.hpp"

#include <algorithm>
#include <vector>

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {
namespace internal_symmetric_linear_multistep_integrator {

using base::make_not_null_unique;

int const startup_step_divisor = 16;

template<typename Position, int order_>
Status SymmetricLinearMultistepIntegrator<Position, order_>::Instance::Solve(
    Instant const& t_final) {
  using Acceleration = typename ODE::Acceleration;
  using Displacement = typename ODE::Displacement;
  using DoubleDisplacement = DoublePrecision<Displacement>;
  using DoubleDisplacements = std::vector<DoubleDisplacement>;
  using DoublePosition = DoublePrecision<Position>;
  using DoublePositions = std::vector<DoublePosition>;

  auto const& ɑ = integrator_.ɑ_;
  auto const& β_numerator = integrator_.β_numerator_;
  auto const& β_denominator = integrator_.β_denominator_;

  if (previous_steps_.size() < order_ - 1) {
    StartupSolve(t_final);
  }

  // Argument checks.
  int const dimension = previous_steps_.back().displacements.size();

  // Time step.
  CHECK_LT(Time(), step_);
  Time const& h = step_;
  // Current time.
  DoublePrecision<Instant> t = previous_steps_.back().time;
  // Order.
  int const k = order_;

  std::vector<Position> positions(dimension);

  DoublePositions Σj_minus_ɑj_qj(dimension);
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
    previous_steps_.pop_front();
    previous_steps_.emplace_back();
    Step& current_step = previous_steps_.back();
    current_step.time = t;
    current_step.accelerations.resize(dimension);

    // Fill the new step.  We skip the division by ɑk as it is equal to 1.0.
    double const ɑk = ɑ[0];
    for (int d = 0; d < dimension; ++d) {
      DoublePosition& current_position = Σj_minus_ɑj_qj[d];
      current_position.Increment(h * h *
                                 Σj_βj_numerator_aj[d] / β_denominator);
      current_step.displacements.push_back(current_position - DoublePosition());
      positions[d] = current_position.value;
      current_state_.positions[d] = current_position;
    }
    equation_.compute_acceleration(t.value,
                                   positions,
                                   current_step.accelerations);

    VelocitySolve(dimension);

    // Inform the caller of the new state.
    current_state_.time = t;
    append_state_(current_state_);
  }

  return Status::OK;
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::
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
    auto* const serialized_step = extension->add_previous_steps();
    for (auto const& displacement : previous_step.displacements) {
      displacement.WriteToMessage(serialized_step->add_displacements());
    }
    for (auto const& acceleration : previous_step.accelerations) {
      acceleration.WriteToMessage(serialized_step->add_accelerations());
    }
    previous_step.time.WriteToMessage(serialized_step->mutable_time());
  }
  integrator_.WriteToMessage(extension->mutable_integrator());
}

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step,
    SymmetricLinearMultistepIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem,
                                             std::move(append_state),
                                             step),
      integrator_(integrator) {
  previous_steps_.emplace_back();
  FillStepFromSystemState(equation_, current_state_, previous_steps_.back());
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::
Instance::StartupSolve(Instant const& t_final) {
  Time const startup_step = step_ / startup_step_divisor;

  CHECK(!previous_steps_.empty());
  CHECK_LT(previous_steps_.size(), order_);

  int startup_step_index = 0;
  auto const startup_append_state =
      [this, &startup_step_index](typename ODE::SystemState const& state) {
        // Stop changing anything once we're done with the startup.  We may be
        // called one more time by the |startup_integrator_|.
        if (previous_steps_.size() < order_) {
          current_state_ = state;
          // The startup integrator has a smaller step.  We do not record all
          // the states it computes, but only those that are a multiple of the
          // main integrator step.
          if (++startup_step_index % startup_step_divisor == 0) {
            CHECK_LT(previous_steps_.size(), order_);
            previous_steps_.emplace_back();
            append_state_(state);
            FillStepFromSystemState(
                equation_, current_state_, previous_steps_.back());
          }
        }
      };

  auto const startup_instance =
      integrator_.startup_integrator_.NewInstance({equation_, &current_state_},
                                                  startup_append_state,
                                                  startup_step);

  startup_instance->Solve(
      std::min(current_state_.time.value +
                   (order_ - previous_steps_.size()) * step_ + step_ / 2.0,
               t_final));

  CHECK_LE(previous_steps_.size(), order_);
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::
Instance::VelocitySolve(int const dimension) {
  using Velocity = typename ODE::Velocity;
  auto const& velocity_integrator = integrator_.velocity_integrator_;

  for (int d = 0; d < dimension; ++d) {
    DoublePrecision<Velocity>& velocity = current_state_.velocities[d];
    auto it = previous_steps_.rbegin();
    Acceleration weighted_acceleration;
    for (int i = 0; i < velocity_integrator.numerators.size; ++i) {
      double const numerator = velocity_integrator.numerators[i];
      weighted_acceleration += numerator * it->accelerations[d];
      ++it;
    }
    velocity.Increment(step_ * weighted_acceleration /
                       velocity_integrator.denominator);
  }
}

template<typename Position, int order_>
void SymmetricLinearMultistepIntegrator<Position, order_>::
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
  equation.compute_acceleration(step.time.value,
                                positions,
                                step.accelerations);
}

template<typename Position, int order_>
SymmetricLinearMultistepIntegrator<Position, order_>::
SymmetricLinearMultistepIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind,
    FixedStepSizeIntegrator<ODE> const& startup_integrator,
    FixedVector<double, half_order_> const & ɑ,
    FixedVector<double, half_order_> const& β_numerator,
    double const β_denominator)
    : FixedStepSizeIntegrator(kind),
      startup_integrator_(startup_integrator),
      velocity_integrator_(AdamsMoultonOrder<velocity_order_>()),
      ɑ_(ɑ),
      β_numerator_(β_numerator),
      β_denominator_(β_denominator) {
  CHECK_EQ(ɑ_[0], 1.0);
  CHECK_EQ(β_numerator_[0], 0.0);
}

template <typename Position, int order_>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
SymmetricLinearMultistepIntegrator<Position, order_>::NewInstance(
    IntegrationProblem<ODE> const& problem,
    typename Integrator<ODE>::AppendState const& append_state,
    Time const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<typename Integrator<ODE>::Instance>(
      new Instance(problem, append_state, step, *this));
}

}  // namespace internal_symmetric_linear_multistep_integrator

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8A() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8A,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -2.0, 2.0, -2.0, 2.0},
      {0.0, 22081.0, -29418.0, 75183.0, -75212.0},
      15120.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const& Quinlan1999Order8B() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_1999_ORDER_8B,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, 0.0, 0.0, -0.5, -1.0},
      {0.0, 192481.0, 6582.0, 816783.0, -156812.0},
      120960.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 8> const&
QuinlanTremaine1990Order8() {
  static SymmetricLinearMultistepIntegrator<Position, 8> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_8,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -2.0, 2.0, -1.0, 0.0},
      {0.0, 17671.0, -23622.0, 61449.0, -50516.0},
      12096.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 10> const&
QuinlanTremaine1990Order10() {
  static SymmetricLinearMultistepIntegrator<Position, 10> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_10,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -1.0, 1.0, -1.0, 1.0, -2.0},
      {0.0, 399187.0, -485156.0, 2391436.0, -2816732.0, 4651330.0},
      241920.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 12> const&
QuinlanTremaine1990Order12() {
  static SymmetricLinearMultistepIntegrator<Position, 12> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_12,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0},
      {0.0,
       90987349.0,
       -229596838.0,
       812627169.0,
       -1628539944.0,
       2714971338.0,
       -3041896548.0},
      53222400.0);
  return integrator;
}

template<typename Position>
SymmetricLinearMultistepIntegrator<Position, 14> const&
QuinlanTremaine1990Order14() {
  static SymmetricLinearMultistepIntegrator<Position, 14> const integrator(
      serialization::FixedStepSizeIntegrator::QUINLAN_TREMAINE_1990_ORDER_14,
      BlanesMoan2002SRKN14A<Position>(),
      {1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0,
       433489274083.0,
       -1364031998256.0,
       5583113380398.0,
       -14154444148720.0,
       28630585332045.0,
       -42056933842656.0,
       48471792742212.0},
      237758976000.0);
  return integrator;
}

}  // namespace integrators
}  // namespace principia
