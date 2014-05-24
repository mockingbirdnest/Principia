
#include <cmath>
#include <ctime>
#include <memory>
#include <vector>

#include "glog/logging.h"

#define TRACE

namespace principia {
namespace integrators {

inline SPRKIntegrator::SPRKIntegrator() : Integrator() {}

inline SPRKIntegrator::~SPRKIntegrator() {}
 
inline std::vector<std::vector<double>> const&
SPRKIntegrator::Order5Optimal() const {
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

inline void SPRKIntegrator::Increment(
      RightHandSideComputation const compute_force,
      AutonomousRightHandSideComputation const compute_velocity,
      Parameters const& parameters,
      Solution* solution) {
  CHECK_NOTNULL(solution);

  std::vector<double> const& a = parameters.coefficients[0];
  //double[] a = coefficients[0];
  std::vector<double> const& b = parameters.coefficients[1];
  //double[] b = coefficients[1];
  int const stages = b.size();
  //int stages = b.Length;
  int const dimension = parameters.q0.size();
  //int dimension = q0.Length;

  // Runge-Kutta time weights.
  std::vector<double> c(stages);
  //double[] c = new double[stages];
  c[0] = 0.0;
  //c[0] = 0;
  // NOTE(phl): Unchanged, possibly suboptimal.
  for (int j = 1; j < stages; ++j) {
    c[j] = c[j - 1] + b[j - 1];
  }
  
  std::vector<double>* p_error;
  std::unique_ptr<std::vector<double>> p_error_deleter;
  if (parameters.p_error == nullptr) {
    p_error = new std::vector<double>(dimension);
  } else {
    p_error = parameters.p_error;
    CHECK_EQ(dimension, p_error->size());
  }
  p_error_deleter.reset(p_error);
  //if (pError == null) {
  //  pError = new double[dimension];
  //}
  
  //TODO(phl): Just say no to code duplication!
  std::vector<double>* q_error;
  std::unique_ptr<std::vector<double>> q_error_deleter;
  if (parameters.q_error == nullptr) {
    q_error = new std::vector<double>(dimension);
  } else {
    q_error = parameters.q_error;
    CHECK_EQ(dimension, q_error->size());
  }
  q_error_deleter.reset(q_error);
  //if (qError == null) {
  //  qError = new double[dimension];
  //}

  double t_error = parameters.t_error;  // NOTE(phl): Don't change |parameters|.

  std::vector<std::vector<double>> Δpstages(stages + 1);
  //double[][] Δpstages = new double[stages + 1][];

  std::vector<std::vector<double>> Δqstages(stages + 1);
  //double[][] Δqstages = new double[stages + 1][];

  for (int i = 0; i < stages + 1; ++i) {
    Δpstages[i].resize(dimension);
    Δqstages[i].resize(dimension);
  }
  //for (int i = 0; i < stages + 1; ++i) {
  //  Δpstages[i] = new double[dimension];
  //  Δqstages[i] = new double[dimension];
  //}

  // Result goes here.
  int const capacity = parameters.sampling_period == 0 ?
    1 : 
    static_cast<int>(
        ceil((((parameters.tmax - parameters.t0) / parameters.Δt) + 1) /
                parameters.sampling_period)) + 1;
  //int capacity = samplingPeriod == 0 ?
  //  1 : (int)Math.Ceiling((((tmax - t0) / Δt) + 1) / samplingPeriod) + 1;
  std::vector<std::vector<double>> q;
  q.reserve(capacity);
  //List<double[]> q = new List<double[]>(capacity);
  std::vector<std::vector<double>> p;
  p.reserve(capacity);
  //List<double[]> p = new List<double[]>(capacity);
  std::vector<double> t;
  t.reserve(capacity);
  //List<double> t = new List<double>(capacity);
  if (parameters.sampling_period != 0) {
    // TODO(phl): Big copies here.
    q.push_back(parameters.q0);
    p.push_back(parameters.p0);
    t.push_back(parameters.t0);
  }
  //if (samplingPeriod != 0) {
  //  q.Add(q0);
  //  p.Add(p0);
  //  t.Add(t0);
  //}

  std::vector<double> q_last(parameters.q0);
  //double[] qLast = (double[])q0.Clone();
  std::vector<double> p_last(parameters.p0);
  //double[] pLast = (double[])p0.Clone();
  double t_last = parameters.t0;
  //double tLast = t0;
  int sampling_phase = 0;
  //int samplingPhase = 0;

  std::vector<double> q_stage(dimension);
  //double[] qStage = new double[dimension];
  std::vector<double> p_stage(dimension);
  //double[] pStage = new double[dimension];
  double tn = parameters.t0;  // Current time.
  //double tn = t0; // Current time.
  double const h = parameters.Δt; // Constant for now.
  //double h = Δt; // Constant for now.
  std::vector<double> f(dimension);  // Current forces.
  //double[] f = new double[dimension]; // Current forces.
  std::vector<double> v(dimension);  // Current velocities.
  //double[] v = new double[dimension]; // Current velocities.

#ifdef TRACE
  int percentage = 0;  // NOTE(phl): Unchanged.
  clock_t running_time = clock();
  //long runningTime = -DateTime.Now.Ticks;
#endif

  // Integration.
    LOG(ERROR)<<tn<<" "<<parameters.tmax;
  while (tn < parameters.tmax) {
  //while (tn < tmax) {
    // Increment SPRK step from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 3.
    for (int k = 0; k < dimension; ++k) {
      Δpstages[0][k] = 0;
      Δqstages[0][k] = 0;
      q_stage[k] = q_last[k];
    }
    //for (int k = 0; k < dimension; ++k) {
    //  Δpstages[0][k] = 0;
    //  Δqstages[0][k] = 0;
    //  qStage[k] = qLast[k];
    //}
    for (int i = 1; i < stages + 1; ++i) {  // NOTE(phl): Unchanged.
      compute_force(tn + c[i - 1] * h, q_stage, &f);
      //computeForce(qStage, tn + c[i - 1] * h, ref f);
      for (int k = 0; k < dimension; ++k) {  // NOTE(phl): Unchanged.
        Δpstages[i][k] = Δpstages[i - 1][k] + h * b[i - 1] * f[k];  // NOTE(phl): Unchanged.
        p_stage[k] = p_last[k] + Δpstages[i][k];
        //pStage[k] = pLast[k] + Δpstages[i][k];
      }
      compute_velocity(p_stage, &v);
      //computeVelocity(pStage, ref v);
      for (int k = 0; k < dimension; ++k) {  // NOTE(phl): Unchanged.
        Δqstages[i][k] = Δqstages[i - 1][k] + h * a[i - 1] * v[k];  // NOTE(phl): Unchanged.
        q_stage[k] = q_last[k] + Δqstages[i][k];
        //qStage[k] = qLast[k] + Δqstages[i][k];
      }
    }
    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {  // NOTE(phl): Unchanged.
      double const Δp = Δpstages[stages][k] + (*p_error)[k];
      //double Δp = Δpstages[stages][k] + pError[k];
      p_stage[k] = p_last[k] + Δp;
      //pStage[k] = pLast[k] + Δp;
      (*p_error)[k] = (p_last[k] - p_stage[k]) + Δp;
      //pError[k] = (pLast[k] - pStage[k]) + Δp;
      p_last[k] = p_stage[k];
      //pLast[k] = pStage[k];
      double Δq = Δqstages[stages][k] + (*q_error)[k];
      //double Δq = Δqstages[stages][k] + qError[k];
      q_stage[k] = q_last[k] + Δq;
      //qStage[k] = qLast[k] + Δq;
      (*q_error)[k] = (q_last[k] - q_stage[k]) + Δq;
      //qError[k] = (qLast[k] - qStage[k]) + Δq;
      q_last[k] = q_stage[k];
      //qLast[k] = qStage[k];
    }

    double const δt = h + t_error;
    //double δt = h + tError;
    tn = tn + δt;  // NOTE(phl): Unchanged.
    t_error = (t_last - tn) + δt;
    //tError = (tLast - tn) + δt;
    t_last = tn;
    //tLast = tn;

    if (parameters.sampling_period != 0) {
    //if (samplingPeriod != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
      //if (samplingPhase % samplingPeriod == 0) {
        t.push_back(tn);
        //t.Add(tn);
        p.push_back(p_stage);
        //p.Add((double[])pStage.Clone());
        q.push_back(q_stage);
        //q.Add((double[])qStage.Clone());
      }
      ++sampling_phase;
      //++samplingPhase;
    }

#ifdef TRACE
    running_time += clock();
    //runningTime += DateTime.Now.Ticks;
    if (floor(tn / parameters.tmax * 100) > percentage) {
    //if (Math.Floor(tn / tmax * 100) > percentage) {
      LOG(ERROR) << "SPRK: " << percentage << "%\ttn = " << tn
                 <<"\tRunning time: " << running_time / (CLOCKS_PER_SEC / 1000)
                 << " ms";
      //Console.WriteLine("SPRK: " + percentage + "%\ttn = " + tn +
      //  "\tRunning time: " + runningTime / 10000 + " ms");
      ++percentage;  // NOTE(phl): Unchanged.
    }
    running_time -= clock();
    //runningTime -= DateTime.Now.Ticks;
#endif
  }
  if (parameters.sampling_period == 0) {
  //if (samplingPeriod == 0) {
    t.push_back(tn);
    //t.Add(tn);
    p.push_back(p_stage);
    //p.Add((double[])pStage.Clone());
    q.push_back(q_stage);
    //q.Add((double[])qStage.Clone());
  }

  // TODO(phl): Maybe avoid all these copies and return ownership?
  for (size_t i = 0; i < p.size(); ++i) {
    solution->momentum[i].quantities = p[i];
    solution->momentum[i].error = (*p_error)[i];
    solution->position[i].quantities = q[i];
    solution->position[i].error = (*q_error)[i];
  }
  solution->time.quantities = t;
  solution->time.error = t_error;

  //return new Solution {
    //momentum = p.ToArray(),
    //momentumError = (double[])pError.Clone(),
    //position = q.ToArray(),
    //positionError = (double[])qError.Clone(),
    //time = t.ToArray(),
    //timeError = tError
  //};
  
}

}  // namespace integrators
}  // namespace principia

#undef TRACE
