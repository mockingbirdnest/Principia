#pragma once

#include "integrators/symplectic_runge_kutta_nystrom_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {

SRKNIntegrator::SRKNIntegrator(std::vector<double> const& a,
                               std::vector<double> const& b)
    : a_(std::move(a)),
      b_(std::move(b)) {
  if (b.front() == 0.0) {
    vanishing_coefficients_ = kFirstBVanishes;
    first_same_as_last_ = std::make_unique<FirstSameAsLast>();
    first_same_as_last_->first = a.front();
    first_same_as_last_->last = a.back();
    a_ = std::vector<double>(a.begin() + 1, a.end());
    b_ = std::vector<double>(b.begin() + 1, b.end());
    a_.back() += first_same_as_last_->first;
    stages_ = b_.size();
    CHECK_EQ(stages_, a_.size());
  } else if (a.back() == 0.0) {
    vanishing_coefficients_ = kLastAVanishes;
    first_same_as_last_ = std::make_unique<FirstSameAsLast>();
    first_same_as_last_->first = b.front();
    first_same_as_last_->last = b.back();
    a_ = std::vector<double>(a.begin(), a.end() - 1);
    b_ = std::vector<double>(b.begin(), b.end() - 1);
    b_.front() += first_same_as_last_->last;
    stages_ = b_.size();
    CHECK_EQ(stages_, a_.size());
  } else {
    vanishing_coefficients_ = kNone;
    first_same_as_last_.reset();
    a_ = a;
    b_ = b;
    stages_ = b_.size();
    CHECK_EQ(stages_, a_.size());
  }

  // Runge-Kutta time weights.
  c_.resize(stages_);
  if (vanishing_coefficients_ == kFirstBVanishes) {
    c_[0] = first_same_as_last_->first;
  } else {
    c_[0] = 0.0;
  }
  for (int j = 1; j < stages_; ++j) {
    c_[j] = c_[j - 1] + a_[j - 1];
  }
}

template<typename Position, typename Velocity,
         typename RightHandSideComputation>
void SRKNIntegrator::SolveTrivialKineticEnergyIncrement(
    RightHandSideComputation compute_acceleration,
    Parameters<Position, Velocity> const& parameters,
    not_null<Solution<Position, Velocity>*> const solution) const {
  switch (vanishing_coefficients_) {
    case kNone:
      SolveTrivialKineticEnergyIncrementOptimized<kNone>(compute_acceleration,
                                                         parameters,
                                                         solution);
      break;
    case kFirstBVanishes:
      SolveTrivialKineticEnergyIncrementOptimized<kFirstBVanishes>(
          compute_acceleration,
          parameters,
          solution);
      break;
    case kLastAVanishes:
      SolveTrivialKineticEnergyIncrementOptimized<kLastAVanishes>(
          compute_acceleration,
          parameters,
          solution);
      break;
    default:
      LOG(FATAL) << "Invalid vanishing coefficients";
  }
}

template<SRKNIntegrator::VanishingCoefficients vanishing_coefficients,
         typename Position, typename Velocity,
         typename RightHandSideComputation>
void SRKNIntegrator::SolveTrivialKineticEnergyIncrementOptimized(
    RightHandSideComputation compute_acceleration,
    SRKNIntegrator::Parameters<Position, Velocity> const& parameters,
    not_null<Solution<Position, Velocity>*> const solution) const {
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
  bool q_and_p_are_synchronized = true;
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

    if (vanishing_coefficients != kNone) {
      should_synchronize = at_end ||
                           (parameters.sampling_period != 0 &&
                            sampling_phase % parameters.sampling_period == 0);
    }

    if (vanishing_coefficients == kFirstBVanishes &&
        q_and_p_are_synchronized) {
      // Desynchronize.
      std::swap(Δqstage_current, Δqstage_previous);
      for (int k = 0; k < dimension; ++k) {
        p_stage[k] = p_last[k].value;
      }
      ADVANCE_ΔQSTAGE(first_same_as_last_->first * h);
      q_and_p_are_synchronized = false;
    }
    for (int i = 0; i < stages_; ++i) {
      std::swap(Δqstage_current, Δqstage_previous);
      std::swap(Δpstage_current, Δpstage_previous);

      // Beware, the p/q order matters here, the two computations depend on one
      // another.

      // By using |tn.error| below we get a time value which is possibly a wee
      // bit more precise.
      if (vanishing_coefficients == kLastAVanishes &&
          q_and_p_are_synchronized && i == 0) {
        ADVANCE_ΔPSTAGE(first_same_as_last_->first * h,
                        tn.value);
        q_and_p_are_synchronized = false;
      } else {
        ADVANCE_ΔPSTAGE(b_[i] * h, tn.value + (tn.error + c_[i] * h));
      }

      if (vanishing_coefficients == kFirstBVanishes &&
          should_synchronize && i == stages_ - 1) {
        ADVANCE_ΔQSTAGE(first_same_as_last_->last * h);
        q_and_p_are_synchronized = true;
      } else {
        ADVANCE_ΔQSTAGE(a_[i] * h);
      }
    }
    if (vanishing_coefficients == kLastAVanishes && should_synchronize) {
      std::swap(Δpstage_current, Δpstage_previous);
      // TODO(egg): the second parameter below is really just tn.value + h.
      ADVANCE_ΔPSTAGE(first_same_as_last_->last * h,
                      tn.value + h);
      q_and_p_are_synchronized = true;
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
