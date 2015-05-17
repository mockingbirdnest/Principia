#pragma once

#include "integrators/implicit_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

#ifdef ADVANCE_ΔQSTAGE
#error ADVANCE_ΔQSTAGE already defined
#else
#define ADVANCE_ΔQSTAGE(step)                                              \
  do {                                                                     \
    Time const step_evaluated = (step);                                    \
    for (int k = 0; k < dimension; ++k) {                                  \
      Displacement const Δq = (*Δqstage_previous)[k] +                     \
                              step_evaluated * v_stage[k];                 \
      q_stage[k] = q_last[k].value + Δq;                                   \
      (*Δqstage_current)[k] = Δq;                                          \
    }                                                                      \
  } while (false)
#endif

#ifdef ADVANCE_ΔVSTAGE
#error ADVANCE_ΔVSTAGE already defined
#else
#define ADVANCE_ΔVSTAGE(step, q_clock)                                     \
  do {                                                                     \
    Time const step_evaluated = (step);                                    \
    compute_acceleration((q_clock), q_stage, &a);                          \
    for (int k = 0; k < dimension; ++k) {                                  \
      Velocity const Δv = (*Δvstage_previous)[k] + step_evaluated * a[k];  \
      v_stage[k] = v_last[k].value + Δv;                                   \
      (*Δvstage_current)[k] = Δv;                                          \
    }                                                                      \
  } while (false)
#endif

namespace principia {

using quantities::Difference;
using quantities::Quotient;

namespace integrators {

template<typename Position>
void IRKNIntegrator::SolveTrivialKineticEnergyIncrement(
    IRKNRightHandSideComputation<Position> compute_acceleration,
    Parameters<Position, Variation<Position>> const& parameters,
    not_null<Solution<Position, Variation<Position>>*> const solution) const {
  using Velocity = Variation<Position>;
  using Displacement = Difference<Position>;
  int const dimension = parameters.initial.positions.size();

  std::vector<std::vector<Displacement>>* Δqstages_current(stages_);
  std::vector<std::vector<Displacement>>* Δqstages_previous(stages_);
  for (int i = 0; i < stages_; ++i) {
    Δqstages_current[i].resize(dimension);
    Δqstages_previous[i].resize(dimension);
  }

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
  std::vector<DoublePrecision<Velocity>> v_last(parameters.initial.momenta);
  int sampling_phase = 0;

  std::vector<Position> q_stage(dimension);
  std::vector<Velocity> v_stage(dimension);
  std::vector<Quotient<Velocity, Time>> a(dimension);  // Current accelerations.

  // The following quantity is generally equal to |Δt|, but during the last
  // iteration, if |tmax_is_exact|, it may differ significantly from |Δt|.
  Time h = parameters.Δt;  // Constant for now.

  // During one iteration of the outer loop below we process the time interval
  // [|tn|, |tn| + |h|[.  |tn| is computed using compensated summation to make
  // sure that we don't have drifts.
  DoublePrecision<Time> tn = parameters.initial.time;

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

    for (int k = 0; k < dimension; ++k) {
      (*Δqstage_current)[k] = Displacement();
      (*Δvstage_current)[k] = Velocity();
      q_stage[k] = q_last[k].value;
    }

    // Fixed-point iteration.
    // NOTE(egg): this does not work for stiff systems, newton iteration is then
    // required.
    for (int i = 0; i < stages_; ++i) {
      Δqstages[i]
    }


    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {
      q_last[k].Increment((*Δqstage_current)[k]);
      v_last[k].Increment((*Δvstage_current)[k]);
      q_stage[k] = q_last[k].value;
      v_stage[k] = v_last[k].value;
    }
    tn.Increment(h);

    if (parameters.sampling_period != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
        solution->emplace_back();
        SystemState<Position, Velocity>* state = &solution->back();
        state->time = tn;
        state->positions.reserve(dimension);
        state->momenta.reserve(dimension);
        for (int k = 0; k < dimension; ++k) {
          state->positions.emplace_back(q_last[k]);
          state->momenta.emplace_back(v_last[k]);
        }
      }
      ++sampling_phase;
    }
  }

  if (parameters.sampling_period == 0) {
    solution->emplace_back();
    SystemState<Position, Velocity>* state = &solution->back();
    state->time = tn;
    state->positions.reserve(dimension);
    state->momenta.reserve(dimension);
    for (int k = 0; k < dimension; ++k) {
      state->positions.emplace_back(q_last[k]);
      state->momenta.emplace_back(v_last[k]);
    }
  }
}

}  // namespace integrators
}  // namespace principia

#undef ADVANCE_ΔQSTAGE
#undef ADVANCE_ΔVSTAGE
