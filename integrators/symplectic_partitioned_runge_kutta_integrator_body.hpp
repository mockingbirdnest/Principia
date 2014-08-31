#pragma once

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "quantities/quantities.hpp"

// Mixed assemblies are not supported by Unity/Mono.
#include "glog/logging.h"

using principia::quantities::Quotient;

namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
inline SPRKIntegrator<Position, Momentum>::SPRKIntegrator() : stages_(0) {}

template<typename Position, typename Momentum>
inline std::vector<std::vector<double>> const&
SPRKIntegrator<Position, Momentum>::Order5Optimal() const {
  static std::vector<std::vector<double>> const order_5_optimal = {
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
       -0.0589796254980311632}};
  return order_5_optimal;
}

template<typename Position, typename Momentum>
inline void SPRKIntegrator<Position, Momentum>::Initialize(
    Coefficients const& coefficients) {
  CHECK_EQ(2, coefficients.size());
  a_ = coefficients[0];
  b_ = coefficients[1];
  stages_ = b_.size();
  CHECK_EQ(stages_, a_.size());

  // Runge-Kutta time weights.
  c_.resize(stages_);
  c_[0] = 0.0;
  for (int j = 1; j < stages_; ++j) {
    c_[j] = c_[j - 1] + b_[j - 1];
  }
}

template<typename Position, typename Momentum>
template<typename AutonomousRightHandSideComputation,
         typename RightHandSideComputation>
void SPRKIntegrator<Position, Momentum>::Solve(
      RightHandSideComputation compute_force,
      AutonomousRightHandSideComputation compute_velocity,
      Parameters const& parameters,
      std::vector<SystemState>* solution) const {
  CHECK_NOTNULL(solution);

  int const dimension = parameters.initial.positions.size();

  std::vector<Position> Δqstage0(dimension);
  std::vector<Position> Δqstage1(dimension);
  std::vector<Momentum> Δpstage0(dimension);
  std::vector<Momentum> Δpstage1(dimension);
  std::vector<Position>* Δqstage_current = &Δqstage1;
  std::vector<Position>* Δqstage_previous = &Δqstage0;
  std::vector<Momentum>* Δpstage_current = &Δpstage1;
  std::vector<Momentum>* Δpstage_previous = &Δpstage0;

  // Dimension the result.
  int const capacity = parameters.sampling_period == 0 ?
    1 :
    static_cast<int>(
        ceil((((parameters.tmax - parameters.initial.time.value) /
                    parameters.Δt) + 1) /
                parameters.sampling_period)) + 1;
  solution->clear();
  solution->reserve(capacity);

  std::vector<DoublePrecision<Position>> q_last(parameters.initial.positions);
  std::vector<DoublePrecision<Momentum>> p_last(parameters.initial.momenta);
  DoublePrecision<Time> t_last = parameters.initial.time;
  int sampling_phase = 0;

  std::vector<Position> q_stage(dimension);
  std::vector<Momentum> p_stage(dimension);
  Time tn = parameters.initial.time.value;  // Current time.
  Time const h = parameters.Δt;  // Constant for now.
  std::vector<Quotient<Momentum, Time>> f(dimension);  // Current forces.
  std::vector<Quotient<Position, Time>> v(dimension);  // Current velocities.

#ifdef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
  int percentage = 0;
  // Initialize |running_time| so that, when we reach the end of the iteration
  // and add clock(), |running_time| will contain the time actually spent in the
  // iteration.
  clock_t running_time = -clock();
#endif

  // Integration.  For details see Wolfram Reference,
  // http://reference.wolfram.com/mathematica/tutorial/NDSolveSPRK.html#74387056
  while (tn < parameters.tmax) {
    // Increment SPRK step from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 3.
    for (int k = 0; k < dimension; ++k) {
      (*Δqstage_current)[k] = Position();
      (*Δpstage_current)[k] = Momentum();
      q_stage[k] = q_last[k].value;
    }
    for (int i = 0; i < stages_; ++i) {
      std::swap(Δqstage_current, Δqstage_previous);
      std::swap(Δpstage_current, Δpstage_previous);
      // Beware, the p/q order matters here, the two computations depend on one
      // another.
      compute_force(tn + c_[i] * h, q_stage, &f);
      for (int k = 0; k < dimension; ++k) {
        Momentum const Δp = (*Δpstage_previous)[k] + h * b_[i] * f[k];
        p_stage[k] = p_last[k].value + Δp;
        (*Δpstage_current)[k] = Δp;
      }
      compute_velocity(p_stage, &v);
      for (int k = 0; k < dimension; ++k) {
        Position const Δq = (*Δqstage_previous)[k] + h * a_[i] * v[k];
        q_stage[k] = q_last[k].value + Δq;
        (*Δqstage_current)[k] = Δq;
      }
    }
    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {
      Position const Δq = (*Δqstage_current)[k] + q_last[k].error;
      q_stage[k] = q_last[k].value + Δq;
      q_last[k].error = (q_last[k].value - q_stage[k]) + Δq;
      q_last[k].value = q_stage[k];
      Momentum const Δp = (*Δpstage_current)[k] + p_last[k].error;
      p_stage[k] = p_last[k].value + Δp;
      p_last[k].error = (p_last[k].value - p_stage[k]) + Δp;
      p_last[k].value = p_stage[k];
    }

    Time const δt = h + t_last.error;
    tn += δt;
    t_last.error = (t_last.value - tn) + δt;
    t_last.value = tn;

    if (parameters.sampling_period != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
        solution->emplace_back();
        SystemState* state = &solution->back();
        state->time = t_last;
        state->positions.reserve(dimension);
        state->momenta.reserve(dimension);
        for (int k = 0; k < dimension; ++k) {
          state->positions.emplace_back(q_last[k]);
          state->momenta.emplace_back(p_last[k]);
        }
      }
      ++sampling_phase;
    }

#ifdef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
    running_time += clock();
    if (floor(tn / parameters.tmax * 100) > percentage) {
      LOG(INFO) << "SPRK: " << percentage << "%\ttn = " << tn
                << "\tRunning time: " << running_time / (CLOCKS_PER_SEC / 1000)
                << " ms";
      ++percentage;
    }
    running_time -= clock();
#endif
  }
  if (parameters.sampling_period == 0) {
    solution->emplace_back();
    SystemState* state = &solution->back();
    state->time = t_last;
    state->positions.reserve(dimension);
    state->momenta.reserve(dimension);
    for (int k = 0; k < dimension; ++k) {
      state->positions.emplace_back(q_last[k]);
      state->momenta.emplace_back(p_last[k]);
    }
  }

#ifdef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
  running_time += clock();
  LOG(INFO) << "SPRK: final running time: "
            << running_time / (CLOCKS_PER_SEC / 1000) << " ms";
#endif
}

}  // namespace integrators
}  // namespace principia
