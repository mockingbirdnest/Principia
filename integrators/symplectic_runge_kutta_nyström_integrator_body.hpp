
#pragma once

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <vector>

#include "geometry/sign.hpp"
#include "numerics/ulp_distance.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using base::make_not_null_unique;
using geometry::Sign;
using numerics::DoublePrecision;
using numerics::ULPDistance;
using quantities::Abs;

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
Status SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                             evaluations, composition>::
Instance::Solve(Instant const& t_final) {
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;
  using Acceleration = typename ODE::Acceleration;

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
  if (integration_direction.Positive()) {
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

  Status status;

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
      status.Update(equation.compute_acceleration(
          t.value + (t.error + c[i] * h), q_stage, g));
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
    append_state(current_state);
  }

  return status;
}

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition> const&
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition>::
Instance::integrator() const {
  return integrator_;
}

template<typename Position,
         int order,
         bool time_reversible,
         int evaluations,
         CompositionMethod composition>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
void SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                           evaluations, composition>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  FixedStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension =
      message
          ->MutableExtension(
              serialization::FixedStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::SymplecticRungeKuttaNystromIntegratorInstance::
                  extension);
}

template<typename Position, int order_, bool time_reversible_, int evaluations_,
         CompositionMethod composition_>
SymplecticRungeKuttaNyströmIntegrator<Position, order_, time_reversible_,
                                      evaluations_, composition_>::
Instance::Instance(IntegrationProblem<ODE> const& problem,
                   AppendState const& append_state,
                   Time const& step,
                   SymplecticRungeKuttaNyströmIntegrator const& integrator)
    : FixedStepSizeIntegrator<ODE>::Instance(problem,
                                             std::move(append_state),
                                             step),
      integrator_(integrator) {}

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition>::
SymplecticRungeKuttaNyströmIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind,
    FixedVector<double, stages_> const& a,
    FixedVector<double, stages_> const& b)
    : FixedStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>>(kind),
      a_(a),
      b_(b) {
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

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition>::NewInstance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(
      new Instance(problem, append_state, step, *this));
}

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition>::
ReadFromMessage(serialization::FixedStepSizeIntegratorInstance const& message,
                IntegrationProblem<ODE> const& problem,
                AppendState const& append_state,
                Time const& step) const {
  CHECK(message.HasExtension(
      serialization::SymplecticRungeKuttaNystromIntegratorInstance::extension))
      << message.DebugString();

  return std::unique_ptr<typename Integrator<ODE>::Instance>(
      new Instance(problem, append_state, step, *this));
}

}  // namespace internal_symplectic_runge_kutta_nyström_integrator

template<typename Method, typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position,
                                      Method::order,
                                      Method::time_reversible,
                                      Method::evaluations,
                                      Method::composition> const&
MakeSymplecticRungeKuttaNyströmIntegrator() {
  static SymplecticRungeKuttaNyströmIntegrator<Position,
                                      Method::order,
                                      Method::time_reversible,
                                      Method::evaluations,
                                      Method::composition> const integrator(
                                      Method::kind, Method::a, Method::b);
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 4, false, 4, BA> const&
McLachlanAtela1992Order4Optimal() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, false, 4, BA> const integrator(
      serialization::FixedStepSizeIntegrator::
          MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL,
      { 0.5153528374311229364,
       -0.085782019412973646,
        0.4415830236164665242,
        0.1288461583653841854},
      { 0.1344961992774310892,
       -0.2248198030794208058,
        0.7563200005156682911,
        0.3340036032863214255});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 4, true, 4, ABA> const&
McLachlan1995SB3A4() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, true, 4, ABA> const integrator(
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SB3A_4,
      { 0.18819521776883821787,
       -0.021528551102171551201,
        0.66666666666666666667,
       -0.021528551102171551201,
        0.18819521776883821787},
      { 0.0,
        1.0,
       -0.5,
       -0.5,
        1.0});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 4, true, 5, ABA> const&
McLachlan1995SB3A5() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, true, 5, ABA> const integrator(
      serialization::FixedStepSizeIntegrator::MCLACHLAN_1995_SB3A_5,
      { 0.4051886183952522772,
       -0.2871440408165240890,
        0.3819554224212718118,
        0.3819554224212718118,
       -0.2871440408165240890,
        0.4051886183952522772},
      { 0.0,
       -0.041095890410958904110,
        0.28813559322033898305,
        0.50592059438123984212,
        0.28813559322033898305,
       -0.041095890410958904110});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 4, true, 6, BAB> const&
BlanesMoan2002SRKN6B() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, true, 6, BAB> const integrator(
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_SRKN_6B,
      { 0.24529895718427100,
        0.60487266571108000,
       -0.35017162289535100,
       -0.35017162289535100,
        0.60487266571108000,
        0.24529895718427100,
        0.0},
      { 0.082984406417405200,
        0.39630980149836800,
       -0.039056304922348600,
        0.1195241940131508,
       -0.039056304922348600,
        0.39630980149836800,
        0.082984406417405200});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 5, false, 6, BA> const&
McLachlanAtela1992Order5Optimal() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 5, false, 6, BA> const integrator(
      serialization::FixedStepSizeIntegrator::
          MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL,
      { 0.339839625839110000,
       -0.088601336903027329,
        0.5858564768259621188,
       -0.603039356536491888,
        0.3235807965546976394,
        0.4423637942197494587},
      { 0.1193900292875672758,
        0.6989273703824752308,
       -0.1713123582716007754,
        0.4012695022513534480,
        0.0107050818482359840,
       -0.0589796254980311632});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 6, true, 7, ABA> const&
OkunborSkeel1994Order6Method13() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 6, true, 7, ABA> const integrator(
      serialization::FixedStepSizeIntegrator::
          OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13,
      {-1.0130879789171747298,
        1.1874295737325427070,
       -0.018335852096460590340,
        0.34399425728109261313,
        0.34399425728109261313,
       -0.018335852096460590340,
        1.1874295737325427070,
       -1.0130879789171747298},
      { 0.0,
        0.00016600692650009894,
       -0.37962421426377360608,
        0.68913741185181063674,
        0.38064159097092574080,
        0.68913741185181063674,
       -0.37962421426377360608,
        0.00016600692650009894});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 6, true, 11, BAB> const&
BlanesMoan2002SRKN11B() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 6, true, 11, BAB> const integrator(
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_SRKN_11B,
      { 0.12322977594627100,
        0.29055379779955800,
       -0.12704921262541700,
       -0.24633176106207500,
        0.35720887279592800,
        0.2047770542914700,
        0.35720887279592800,
       -0.24633176106207500,
       -0.12704921262541700,
        0.29055379779955800,
        0.12322977594627100,
        0.0},
      { 0.041464998518262400,
        0.19812867191806700,
       -0.040006192104153300,
        0.075253984301580700,
       -0.011511387420687900,
        0.23666992478693110,
        0.23666992478693110,
       -0.011511387420687900,
        0.075253984301580700,
       -0.040006192104153300,
        0.19812867191806700,
        0.041464998518262400});
  return integrator;
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 6, true, 14, ABA> const&
BlanesMoan2002SRKN14A() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 6, true, 14, ABA> const integrator(
      serialization::FixedStepSizeIntegrator::BLANES_MOAN_2002_SRKN_14A,
      { 0.037859319840611600,
        0.10263563310243500,
       -0.025867888266558700,
        0.31424140307144700,
       -0.13014445951741500,
        0.10641770036954300,
       -0.0087942431285105800,
        0.2073050690568954,
       -0.0087942431285105800,
        0.10641770036954300,
       -0.13014445951741500,
        0.31424140307144700,
       -0.025867888266558700,
        0.10263563310243500,
        0.037859319840611600},
      { 0.0,
        0.091719152624461650,
        0.18398317000500600,
       -0.056534365832888270,
        0.0049146887747128540,
        0.14376112716835800,
        0.32856769374680400,
       -0.19641146648645423,
       -0.19641146648645423,
        0.32856769374680400,
        0.14376112716835800,
        0.0049146887747128540,
       -0.056534365832888270,
        0.18398317000500600,
        0.091719152624461650});
  return integrator;
}

}  // namespace integrators
}  // namespace principia
