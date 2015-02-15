#pragma once

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "quantities/quantities.hpp"

// Mixed assemblies are not supported by Unity/Mono.
#include "glog/logging.h"

namespace principia {

using quantities::Quotient;

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
      not_null<std::vector<SystemState>*> const solution) const {
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
  int sampling_phase = 0;

  std::vector<Position> q_stage(dimension);
  std::vector<Momentum> p_stage(dimension);
  std::vector<Quotient<Momentum, Time>> f(dimension);  // Current forces.
  std::vector<Quotient<Position, Time>> v(dimension);  // Current velocities.

  // The following quantity is generally equal to |Δt|, but during the last
  // iteration, if |tmax_is_exact|, it may differ significantly from |Δt|.
  Time h = parameters.Δt;  // Constant for now.

  // During one iteration of the outer loop below we process the time interval
  // [|tn|, |tn| + |h|[.  |tn| is computed using compensated summation to make
  // sure that we don't have drifts.
  DoublePrecision<Time> tn = parameters.initial.time;

#ifdef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
  int percentage = 0;
  // Initialize |running_time| so that, when we reach the end of the iteration
  // and add clock(), |running_time| will contain the time actually spent in the
  // iteration.
  clock_t running_time = -clock();
#endif

  // Integration.  For details see Wolfram Reference,
  // http://reference.wolfram.com/mathematica/tutorial/NDSolveSPRK.html#74387056
  bool at_end = !parameters.tmax_is_exact && parameters.tmax < tn.value + h;
  while (!at_end) {
    // Check if this is the last interval and if so process it appropriately.
    if (parameters.tmax_is_exact) {
      // If |tn| is getting close to |tmax|, use |tmax| as the upper bound of
      // the interval and update |h| accordingly.  The bound chosen here for
      // |tmax| ensures that we don't end up with a ridiculously small last
      // interval: we'd rather make the last interval a bit bigger.  More
      // precisely, the last interval generally has a length between 0.5 Δt and
      // 1.5 Δt, unless it is also the first interval.
      // NOTE(phl): This may lead to convergence as bad as (1.5 Δt)^5 rather
      // than Δt^5.
      if (parameters.tmax <= tn.value + 3 * h / 2) {
        at_end = true;
        h = (parameters.tmax - tn.value) - tn.error;
      }
    } else if (parameters.tmax < tn.value + 2 * h) {
      // If the next interval would overshoot, make this the last interval but
      // stick to the same step.
      at_end = true;
    }
    // Here |h| is the length of the current time interval and |tn| is its
    // start.

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

      // By using |tn.error| below we get a time value which is possibly a wee
      // bit more precise.
      compute_force(tn.value + (tn.error + c_[i] * h), q_stage, &f);

      // Beware, the p/q order matters here, the two computations depend on one
      // another.
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
      q_last[k].Increment((*Δqstage_current)[k]);
      p_last[k].Increment((*Δpstage_current)[k]);
      q_stage[k] = q_last[k].value;
      p_stage[k] = p_last[k].value;
    }
    tn.Increment(h);

    if (parameters.sampling_period != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
        solution->emplace_back();
        SystemState* state = &solution->back();
        state->time = tn;
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
    while (floor(100 * (tn.value - parameters.initial.time.value) /
                 (parameters.tmax - parameters.initial.time.value)) >
           percentage) {
      LOG(INFO) << "SPRK: " << percentage << "%\ttn = " << tn.value
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
    state->time = tn;
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
