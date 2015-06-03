#pragma once

#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

#include <vector>

#include "geometry/sign.hpp"

namespace principia {

using geometry::Sign;

namespace integrators {

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                      evaluations, composition>::
SymplecticRungeKuttaNyströmIntegrator(FixedVector<double, stages_> const& a,
                                      FixedVector<double, stages_> const& b)
    : a_(a),
      b_(b) {
  DoublePrecision<double> c_i = 0.0;
  for (int i = 0; i < stages_; ++i) {
    c_[i] = c_i.value;
    c_i.Increment(a_[i]);
  }
  CHECK_EQ(1.0, c_i.value);
  if (composition == kABA) {
    CHECK_EQ(0.0, b_[0]);
  } else if (composition == kBAB) {
    CHECK_EQ(0.0, a_[stages_ - 1]);
  }
  if (time_reversible) {
    switch (composition) {
      case kABA:
        for (int i = 0; i < stages_; ++i) {
          CHECK_EQ(a_[i], a_[stages_ - 1 - i]);
        }
        for (int i = 0; i < stages_ - 1; ++i) {
          CHECK_EQ(b_[i + 1], b_[stages_ - 1 - i]);
        }
        break;
       case kBAB:
        for (int i = 0; i < stages_ - 1; ++i) {
          CHECK_EQ(b_[i], a_[stages_ - 2 - i]);
        }
        for (int i = 0; i < stages_; ++i) {
          CHECK_EQ(b_[i], b_[stages_ - 1 - i]);
        }
        break;
      case kBA:
        LOG(FATAL) << "Time-reversible compositions have the FSAL property";
        break;
      default:
        LOG(FATAL) << "Invalid CompositionMethod";
    }
  }
}

template<typename Position, int order, bool time_reversible, int evaluations,
         CompositionMethod composition>
void SymplecticRungeKuttaNyströmIntegrator<Position, order, time_reversible,
                                           evaluations, composition>::Solve(
    IntegrationProblem<ODE> const& problem,
    Time const& step) const {
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;
  using Acceleration = typename ODE::Acceleration;

  // Argument checks.
  CHECK_NOTNULL(problem.initial_state);
  int const dimension = problem.initial_state->positions.size();
  CHECK_EQ(dimension, problem.initial_state->velocities.size());
  CHECK_NE(Time(), step);
  Sign const integration_direction = Sign(step);
  if (integration_direction.Positive()) {
    // Integrating forward.
    CHECK_LT(problem.initial_state->time.value, problem.t_final);
  } else {
    // Integrating backward.
    CHECK_GT(problem.initial_state->time.value, problem.t_final);
  }

  typename ODE::SystemState current_state = *problem.initial_state;

  // Time step.
  Time h = step;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment.
  std::vector<Displacement> ∆q(dimension);
  // Velocity increment.
  std::vector<Velocity> ∆v(dimension);
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

  bool at_end = false;

  for (;;) {
    Time const time_to_end = (problem.t_final - t.value) - t.error;
    at_end = integration_direction * h > integration_direction * time_to_end;
    if (at_end) {
      break;
    }

    std::fill(∆q.begin(), ∆q.end(), Displacement{});
    std::fill(∆v.begin(), ∆v.end(), Velocity{});
    for (int i = 0; i < stages_; ++i) {
      for (int k = 0; k < dimension; ++k) {
        q_stage[k] = q[k].value + ∆q[k];
      }
      problem.equation.compute_acceleration(t.value + c_[i] * h, q_stage, &g);
      for (int k = 0; k < dimension; ++k) {
        // TODO(egg): reformulate to reduce roundoff error.
        ∆v[k] += h * b_[i] * g[k];
        ∆q[k] += h * a_[i] * (v[k].value + ∆v[k]);
      }
    }
    // TODO(egg): handle the FSAL case (types ABA and BAB).

    // Increment the solution.
    t.Increment(h);
    for (int k = 0; k < dimension; ++k) {
      q[k].Increment(∆q[k]);
      v[k].Increment(∆v[k]);
    }
    problem.append_state(current_state);
  }
}

template<typename Position>
SymplecticRungeKuttaNyströmIntegrator<Position, 4, false, 4, kBA> const&
McLachlanAtela1992Order4Optimal() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, false, 4, kBA> const integrator(
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
SymplecticRungeKuttaNyströmIntegrator<Position, 4, true, 4, kABA> const&
McLachlan1995SB3A4() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, true, 4, kABA> const integrator(
      { 0.18819521776883821787,
        0.021528551102171551201,
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
SymplecticRungeKuttaNyströmIntegrator<Position, 4, true, 5, kABA> const&
McLachlan1995SB3A5() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, true, 5, kABA> const integrator(
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
SymplecticRungeKuttaNyströmIntegrator<Position, 4, true, 6, kBAB> const&
BlanesMoan2002SRKN6B() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 4, true, 6, kBAB> const integrator(
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
SymplecticRungeKuttaNyströmIntegrator<Position, 5, false, 6, kBA> const&
McLachlanAtela1992Order5Optimal() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 5, false, 6, kBA> const integrator(
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
SymplecticRungeKuttaNyströmIntegrator<Position, 6, true, 7, kABA> const&
OkunborSkeel1994Order6Method13() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 6, true, 7, kABA> const integrator(
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
SymplecticRungeKuttaNyströmIntegrator<Position, 6, true, 11, kBAB> const&
BlanesMoan2002SRKN11B() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 6, true, 11, kBAB> const integrator(
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
SymplecticRungeKuttaNyströmIntegrator<Position, 6, true, 14, kABA> const&
BlanesMoan2002SRKN14A() {
  static SymplecticRungeKuttaNyströmIntegrator<
             Position, 6, true, 14, kABA> const integrator(
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
