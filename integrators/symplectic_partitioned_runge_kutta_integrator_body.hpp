
#include <algorithm>
#include <cmath>
#include <ctime>
#include <memory>
#include <vector>

#include "quantities/quantities.hpp"

// Mixed assemblies are not supported by Unity/Mono.
#ifndef _MANAGED
#include "glog/logging.h"
#endif

using principia::quantities::Quotient;

namespace principia {
namespace integrators {

namespace {

template <typename T>
inline std::vector<T>* PointerOrNew(int const dimension,
                                    std::vector<T>* const in) {
  if (in == nullptr) {
    return new std::vector<T>(dimension);
  } else {
#ifndef _MANAGED
    CHECK_EQ(dimension, in->size());
#endif
    return in;
  }
}

}  // namespace

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
#ifndef _MANAGED
  CHECK_EQ(2, coefficients.size());
#endif
  a_ = coefficients[0];
  b_ = coefficients[1];
  stages_ = b_.size();
#ifndef _MANAGED
  CHECK_EQ(stages_, a_.size());
#endif

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
      Solution* solution) const {
#ifndef _MANAGED
  CHECK_NOTNULL(solution);
#endif

  int const dimension = parameters.q0.size();
  std::unique_ptr<std::vector<Position>> q_error(
      PointerOrNew(dimension, parameters.q_error));
  std::unique_ptr<std::vector<Momentum>> p_error(
      PointerOrNew(dimension, parameters.p_error));
  Time t_error = parameters.t_error;

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
        ceil((((parameters.tmax - parameters.t0) / parameters.Δt) + 1) /
                parameters.sampling_period)) + 1;
  solution->time.quantities.clear();
  solution->time.quantities.reserve(capacity);
  solution->momentum.resize(dimension);
  solution->position.resize(dimension);
  for (int k = 0; k < dimension; ++k) {
    solution->position[k].quantities.clear();
    solution->position[k].quantities.reserve(capacity);
    solution->momentum[k].quantities.clear();
    solution->momentum[k].quantities.reserve(capacity);
  }

  std::vector<Position> q_last(parameters.q0);
  std::vector<Momentum> p_last(parameters.p0);
  Time t_last = parameters.t0;
  int sampling_phase = 0;

  std::vector<Position> q_stage(dimension);
  std::vector<Momentum> p_stage(dimension);
  Time tn = parameters.t0;  // Current time.
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
      (*Δqstage_current)[k] = Length();
      (*Δpstage_current)[k] = Momentum();
      q_stage[k] = q_last[k];
    }
    for (int i = 0; i < stages_; ++i) {
      std::swap(Δqstage_current, Δqstage_previous);
      std::swap(Δpstage_current, Δpstage_previous);
      // Beware, the p/q order matters here, the two computations depend on one
      // another.
      compute_force(tn + c_[i] * h, q_stage, &f);
      for (int k = 0; k < dimension; ++k) {
        Momentum const Δp = (*Δpstage_previous)[k] + h * b_[i] * f[k];
        p_stage[k] = p_last[k] + Δp;
        (*Δpstage_current)[k] = Δp;
      }
      compute_velocity(p_stage, &v);
      for (int k = 0; k < dimension; ++k) {
        Position const Δq = (*Δqstage_previous)[k] + h * a_[i] * v[k];
        q_stage[k] = q_last[k] + Δq;
        (*Δqstage_current)[k] = Δq;
      }
    }
    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {
      Position const Δq = (*Δqstage_current)[k] + (*q_error)[k];
      q_stage[k] = q_last[k] + Δq;
      (*q_error)[k] = (q_last[k] - q_stage[k]) + Δq;
      q_last[k] = q_stage[k];
      Momentum const Δp = (*Δpstage_current)[k] + (*p_error)[k];
      p_stage[k] = p_last[k] + Δp;
      (*p_error)[k] = (p_last[k] - p_stage[k]) + Δp;
      p_last[k] = p_stage[k];
    }

    Time const δt = h + t_error;
    tn += δt;
    t_error = (t_last - tn) + δt;
    t_last = tn;

    if (parameters.sampling_period != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
        solution->time.quantities.push_back(tn);
        for (int k = 0; k < dimension; ++k) {
          solution->position[k].quantities.push_back(q_stage[k]);
          solution->momentum[k].quantities.push_back(p_stage[k]);
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
    solution->time.quantities.push_back(tn);
    for (int k = 0; k < dimension; ++k) {
      solution->position[k].quantities.push_back(q_stage[k]);
      solution->momentum[k].quantities.push_back(p_stage[k]);
    }
  }

  solution->time.error = t_error;
  for (int k = 0; k < dimension; ++k) {
    solution->position[k].error = (*q_error)[k];
    solution->momentum[k].error = (*p_error)[k];
  }

#ifdef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
  running_time += clock();
  LOG(INFO) << "SPRK: final running time: "
            << running_time / (CLOCKS_PER_SEC / 1000) << " ms";
#endif
}

}  // namespace integrators
}  // namespace principia
