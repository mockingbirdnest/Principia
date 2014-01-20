using Integrators;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NewtonianPhysics {
  public class NBodySystem {
    public const double G = 6.67384E-11;
    public Body[] bodies;
    public NBodySystem(Body[] bodies, double t0, double tError = 0) {
      this.bodies = bodies;
      dimension = bodies.Length * 3;
      q = new double[dimension];
      v = new double[dimension];
      qError = new double[dimension];
      vError = new double[dimension];
      t = t0;
      this.tError = tError;
      updateState();
    }
    public void Evolve(double tmax, double maxTimestep) {
      SymplecticPartitionedRungeKutta.Solution solution
        = Integrators.SymplecticPartitionedRungeKutta.IncrementSPRK(
            computeForce: computeAccelerations,
            computeVelocity: computeVelocities,
            q0: q, p0: v, t0: t, tmax: tmax,
            Δt: (tmax - t) / Math.Ceiling((tmax - t) / maxTimestep),
            coefficients: SymplecticPartitionedRungeKutta.Order5Optimal,
            samplingPeriod: 0, qError: qError, pError: vError, tError: tError);
      q = solution.position[0];
      v = solution.momentum[0];
      t = solution.time[0];
      qError = solution.positionError;
      vError = solution.momentumError;
      tError = solution.timeError;
      updateBodies();
    }
    public void RecalculateAllPredictions(double tmax,
                                         double maxTimestep,
                                         int samplingPeriod) {
      SymplecticPartitionedRungeKutta.Solution solution
        = Integrators.SymplecticPartitionedRungeKutta.IncrementSPRK(
            computeForce: computeAccelerations,
            computeVelocity: computeVelocities,
            q0: q, p0: v, t0: t, tmax: tmax,
            Δt: (tmax - t) / Math.Ceiling((tmax - t) / maxTimestep),
            coefficients: SymplecticPartitionedRungeKutta.Order5Optimal,
            samplingPeriod: samplingPeriod, qError: qError, pError: vError,
            tError: tError);
      for (int b = 0; b < bodies.Length; ++b) {
        bodies[b].predictedTrajectory.Clear();
      }
      for (int i = 0; i < solution.position.Length; ++i) {
        for (int b = 0; b < bodies.Length; ++b) {
          bodies[b].predictedTrajectory.Add(new Event {
            q = new SpatialCoordinates {
              x = solution.position[i][3 * b],
              y = solution.position[i][3 * b + 1],
              z = solution.position[i][3 * b + 2]
            },
            t = solution.time[i]
          });
        }
      }
    }

    #region private

    private int dimension;
    private double[] q;
    private double[] qError;
    private double[] qPredicted;
    private double[] qPredictedError;
    private double t;
    private double tError;
    private double[] v;
    private double[] vError;
    private double[] vPredicted;
    private double[] vPredictedError;
    private void computeAccelerations(double[] q,
                                      double t,
                                      ref double[] result) {
      for (int k = 0; k < dimension; ++k) {
        result[k] = 0;
      }

      for (int b1 = 0; b1 < bodies.Length; ++b1) {
        for (int b2 = b1 + 1; b2 < bodies.Length; ++b2) {
          if (!(bodies[b1].massless && bodies[b2].massless)) {
            double Δq0 = q[3 * b1] - q[3 * b2];
            double Δq1 = q[3 * b1 + 1] - q[3 * b2 + 1];
            double Δq2 = q[3 * b1 + 2] - q[3 * b2 + 2];

            double squaredDistance = Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
            double denominator = squaredDistance * Math.Sqrt(squaredDistance);

            if (!bodies[b2].massless) {
              double μ2OverRSquared
                = bodies[b2].gravitationalParameter / denominator;
              result[3 * b1] -= Δq0 * μ2OverRSquared;
              result[3 * b1 + 1] -= Δq1 * μ2OverRSquared;
              result[3 * b1 + 2] -= Δq2 * μ2OverRSquared;
            }
            // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
            // sive corporum duorum actiones in se mutuo semper esse æquales &
            // in partes contrarias dirigi.
            if (!bodies[b1].massless) {
              double μ1OverRSquared
                = bodies[b1].gravitationalParameter / denominator;
              result[3 * b2] += Δq0 * μ1OverRSquared;
              result[3 * b2 + 1] += Δq1 * μ1OverRSquared;
              result[3 * b2 + 2] += Δq2 * μ1OverRSquared;
            }
          }
        }
      }
    }
    private void computeVelocities(double[] v, ref double[] result) {
      result = v;
    }
    private void updateBodies() {
      for (int b = 0; b < bodies.Length; ++b) {
        bodies[b].q.x = q[3 * b];
        bodies[b].q.y = q[3 * b + 1];
        bodies[b].q.z = q[3 * b + 2];
        bodies[b].v.x = v[3 * b];
        bodies[b].v.y = v[3 * b + 1];
        bodies[b].v.z = v[3 * b + 2];
        bodies[b].qError.x = qError[3 * b];
        bodies[b].qError.y = qError[3 * b + 1];
        bodies[b].qError.z = qError[3 * b + 2];
        bodies[b].vError.x = vError[3 * b];
        bodies[b].vError.y = vError[3 * b + 1];
        bodies[b].vError.z = vError[3 * b + 2];
      }
    }
    private void updateState() {
      for (int b = 0; b < bodies.Length; ++b) {
        q[3 * b] = bodies[b].q.x;
        q[3 * b + 1] = bodies[b].q.y;
        q[3 * b + 2] = bodies[b].q.z;
        v[3 * b] = bodies[b].v.x;
        v[3 * b + 1] = bodies[b].v.y;
        v[3 * b + 2] = bodies[b].v.z;
        qError[3 * b] = bodies[b].qError.x;
        qError[3 * b + 1] = bodies[b].qError.y;
        qError[3 * b + 2] = bodies[b].qError.z;
        vError[3 * b] = bodies[b].vError.x;
        vError[3 * b + 1] = bodies[b].vError.y;
        vError[3 * b + 2] = bodies[b].vError.z;
      }
    }

    #endregion private
  }
}