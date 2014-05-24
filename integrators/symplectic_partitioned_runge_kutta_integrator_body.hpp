
#include <cmath>
#include <memory>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace integrators {

SPRKIntegrator::SPRKIntegrator() : Integrator() {}

SPRKIntegrator::~SPRKIntegrator() {}

void SPRKIntegrator::Increment(
      RightHandSideComputation const compute_force,
      AutonomousRightHandSideComputation const compute_velocity,
      Parameters const& parameters,
      Solution* solution) {
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
  //----
  double tn = t0; // Current time.
  double h = Δt; // Constant for now.
  double[] f = new double[dimension]; // Current forces.
  double[] v = new double[dimension]; // Current velocities.

#if TRACE
  int percentage = 0;
  long runningTime = -DateTime.Now.Ticks;
#endif

  // Integration.
  while (tn < tmax) {
    // Increment SPRK step from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 3.
    for (int k = 0; k < dimension; ++k) {
      Δpstages[0][k] = 0;
      Δqstages[0][k] = 0;
      qStage[k] = qLast[k];
    }
    for (int i = 1; i < stages + 1; ++i) {
      computeForce(qStage, tn + c[i - 1] * h, ref f);
      for (int k = 0; k < dimension; ++k) {
        Δpstages[i][k] = Δpstages[i - 1][k] + h * b[i - 1] * f[k];
        pStage[k] = pLast[k] + Δpstages[i][k];
      }
      computeVelocity(pStage, ref v);
      for (int k = 0; k < dimension; ++k) {
        Δqstages[i][k] = Δqstages[i - 1][k] + h * a[i - 1] * v[k];
        qStage[k] = qLast[k] + Δqstages[i][k];
      }
    }
    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {
      double Δp = Δpstages[stages][k] + pError[k];
      pStage[k] = pLast[k] + Δp;
      pError[k] = (pLast[k] - pStage[k]) + Δp;
      pLast[k] = pStage[k];
      double Δq = Δqstages[stages][k] + qError[k];
      qStage[k] = qLast[k] + Δq;
      qError[k] = (qLast[k] - qStage[k]) + Δq;
      qLast[k] = qStage[k];
    }

    double δt = h + tError;
    tn = tn + δt;
    tError = (tLast - tn) + δt;
    tLast = tn;

    if (samplingPeriod != 0) {
      if (samplingPhase % samplingPeriod == 0) {
        t.Add(tn);
        p.Add((double[])pStage.Clone());
        q.Add((double[])qStage.Clone());
      }
      ++samplingPhase;
    }

#if TRACE
    runningTime += DateTime.Now.Ticks;
    if (Math.Floor(tn / tmax * 100) > percentage) {
      Console.WriteLine("SPRK: " + percentage + "%\ttn = " + tn +
        "\tRunning time: " + runningTime / 10000 + " ms");
      ++percentage;
    }
    runningTime -= DateTime.Now.Ticks;
#endif
  }
  if (samplingPeriod == 0) {
    t.Add(tn);
    p.Add((double[])pStage.Clone());
    q.Add((double[])qStage.Clone());
  }
  return new Solution {
    momentum = p.ToArray(),
    momentumError = (double[])pError.Clone(),
    position = q.ToArray(),
    positionError = (double[])qError.Clone(),
    time = t.ToArray(),
    timeError = tError
  };
}

}  // namespace integrators
}  // namespace principia
