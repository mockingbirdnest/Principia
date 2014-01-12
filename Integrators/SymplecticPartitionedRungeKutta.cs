using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Integrators {
  public static class SymplecticPartitionedRungeKutta {
    // From "The accuracy of symplectic integrators", Robert I. McLachlan and
    // Pau Atela (1992).
    public static double[][] Order5Optimal = {
      new double[] { 0.339839625839110000,
                    -0.088601336903027329,
                     0.5858564768259621188,
                     -0.603039356536491888,
                     0.3235807965546976394,
                     0.4423637942197494587 },
      new double[] { 0.1193900292875672758,
                     0.6989273703824752308,
                    -0.1713123582716007754,
                     0.4012695022513534480,
                     0.0107050818482359840,
                    -0.0589796254980311632 } };
    public static Solution IncrementSPRK(ComputeRightHandSide computeForce,
                      ComputeAutonomousRightHandSide computeVelocity,
                      double[] q0, double[] p0,
                      double t0, double tmax,
                      double Δt,
                      double[][] coefficients,
                      int samplingPeriod,
                      double[] qError = null,
                      double[] pError = null,
                      double tError = 0) {
      double[] a = coefficients[0];
      double[] b = coefficients[1];
      int stages = b.Length;
      int dimension = q0.Length;
      // Runge-Kutta time weights.
      double[] c = new double[stages];
      c[0] = 0;
      for (int j = 1; j < stages; ++j) {
        c[j] = c[j - 1] + b[j - 1];
      }

      if (pError == null) {
        pError = new double[dimension];
      }
      if (qError == null) {
        qError = new double[dimension];
      }

      double[][] Δpstages = new double[stages + 1][];
      double[][] Δqstages = new double[stages + 1][];
      for (int i = 0; i < stages + 1; ++i) {
        Δpstages[i] = new double[dimension];
        Δqstages[i] = new double[dimension];
      }

      // Result goes here.
      List<double[]> q = new List<double[]>(
          (int)Math.Ceiling((((tmax - t0) / Δt) + 1) / samplingPeriod) + 1
        );
      List<double[]> p = new List<double[]>(
          (int)Math.Ceiling((((tmax - t0) / Δt) + 1) / samplingPeriod) + 1
        );
      List<double> t = new List<double>(
          (int)Math.Ceiling((((tmax - t0) / Δt) + 1) / samplingPeriod) + 1
        );
      q.Add(q0);
      p.Add(p0);
      t.Add(t0);

      double[] qLast = (double[])q0.Clone();
      double[] pLast = (double[])p0.Clone();
      double tLast = t0;
      int samplingPhase = 0;

      double[] qStage = new double[dimension];
      double[] pStage = new double[dimension];
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
        // TODO: choose timestep here.

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
    public struct Solution {
      public double[][] momentum;
      public double[] momentumError;
      public double[][] position;
      public double[] positionError;
      public double[] time;
      public double timeError;
    }
  }
}