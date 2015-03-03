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
inline typename SPRKIntegrator<Position, Momentum>::Coefficients const&
SPRKIntegrator<Position, Momentum>::Leapfrog() const {
  static Coefficients const leapfrog = {{ 0.5, 0.5}, { 0.0, 1.0}};
  return leapfrog;
}

template<typename Position, typename Momentum>
inline typename SPRKIntegrator<Position, Momentum>::Coefficients const&
SPRKIntegrator<Position, Momentum>::Order4FirstSameAsLast() const {
  static Coefficients const order_4_first_same_as_last = {
      { 0.6756035959798288170,
       -0.1756035959798288170,
       -0.1756035959798288170,
        0.6756035959798288170},
      { 0.0,
        1.351207191959657634,
       -1.702414383919315268,
        1.351207191959657634}};
  return order_4_first_same_as_last;
}

template<typename Position, typename Momentum>
inline typename SPRKIntegrator<Position, Momentum>::Coefficients const&
SPRKIntegrator<Position, Momentum>::Order5Optimal() const {
  static Coefficients const order_5_optimal = {
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
  if (coefficients[1].front() == 0.0) {
    vanishing_coefficients_ = FirstBVanishes;
    first_same_as_last_ = std::make_unique<FirstSameAsLast>();
    first_same_as_last_->first = coefficients[0].front();
    first_same_as_last_->last = coefficients[0].back();
    a_ = std::vector<double>(coefficients[0].begin() + 1,
                             coefficients[0].end());
    b_ = std::vector<double>(coefficients[0].begin() + 1,
                             coefficients[0].end());
    a_.back() += first_same_as_last_->first;
    stages_ = b_.size();
    CHECK_EQ(stages_, a_.size());
  } else if (coefficients[0].back() == 0.0) {
    vanishing_coefficients_ = LastAVanishes;
    first_same_as_last_ = std::make_unique<FirstSameAsLast>();
    first_same_as_last_->first = coefficients[1].front();
    first_same_as_last_->last = coefficients[1].back();
    a_ = std::vector<double>(coefficients[0].begin(),
                             coefficients[0].end() - 1);
    b_ = std::vector<double>(coefficients[0].begin(),
                             coefficients[0].end() - 1);
    b_.front() += first_same_as_last_->last;
    stages_ = b_.size();
    CHECK_EQ(stages_, a_.size());
  } else {
    a_ = coefficients[0];
    b_ = coefficients[1];
    stages_ = b_.size();
    CHECK_EQ(stages_, a_.size());
  }

  // Runge-Kutta time weights.
  c_.resize(stages_);
  if (vanishing_coefficients_ == FirstBVanishes) {
    c_[0] = first_same_as_last_->first;
  } else {
    c_[0] = 0.0;
  }
  for (int j = 1; j < stages_; ++j) {
    c_[j] = c_[j - 1] + a_[j - 1];
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
  switch (vanishing_coefficients_) {
    case None:
      SolveOptimized<None>(
          compute_force, compute_velocity, parameters, solution);
      break;
    case FirstBVanishes:
      SolveOptimized<FirstBVanishes>(
          compute_force, compute_velocity, parameters, solution);
      break;
    case LastAVanishes:
      SolveOptimized<LastAVanishes>(
          compute_force, compute_velocity, parameters, solution);
      break;
    default:
      LOG(FATAL) << "Invalid vanishing coefficients";
  }
}

template<typename Position, typename Momentum>
template<VanishingCoefficients vanishing_coefficients,
         typename AutonomousRightHandSideComputation,
         typename RightHandSideComputation>
void SPRKIntegrator<Position, Momentum>::SolveOptimized(
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

  // Whether position and momentum are synchronized between steps, relevant for
  // first-same-as-last (FSAL) integrators. Time is always synchronous with
  // position.
  bool synchronized = true;
  bool should_synchronize = false;

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

    if (first_same_as_last) {
      if (synchronized) {
        // Desynchronize.
        for (int k = 0; k < dimension; ++k) {
          p_stage[k] = p_last[k].value;
        }
        compute_velocity(p_stage, &v);
        for (int k = 0; k < dimension; ++k) {
          Position const Δq = (*Δqstage_previous)[k] +
                              h * first_same_as_last_->a_first * v[k];
          q_stage[k] = q_last[k].value + Δq;
          (*Δqstage_current)[k] = Δq;
        }
      }
      should_synchronize = at_end ||
                           (parameters.sampling_period != 0 &&
                            sampling_phase % parameters.sampling_period == 0);
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
      if (first_same_as_last && should_synchronize && i == stages_ - 1) {
        break;
      }
      for (int k = 0; k < dimension; ++k) {
        Position const Δq = (*Δqstage_previous)[k] + h * a_[i] * v[k];
        q_stage[k] = q_last[k].value + Δq;
        (*Δqstage_current)[k] = Δq;
      }
    }
    if (first_same_as_last && should_synchronize) {
      for (int k = 0; k < dimension; ++k) {
        Position const Δq = (*Δqstage_previous)[k] +
                            h * first_same_as_last_->a_last * v[k];
        q_stage[k] = q_last[k].value + Δq;
        (*Δqstage_current)[k] = Δq;
      }
      synchronized = true;
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

}

}  // namespace integrators
}  // namespace principia
