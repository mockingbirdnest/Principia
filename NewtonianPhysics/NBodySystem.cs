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
    public NBodySystem(Body[] bodies) {
      this.bodies = bodies;
      dimension = bodies.Length * 3;
      q = new double[dimension];
      v = new double[dimension];
      qError = new double[dimension];
      vError = new double[dimension];
      updateState();
    }
    public void Evolve(double duration, double timestep) {
      SymplecticPartitionedRungeKutta.Solution solution
        = Integrators.SymplecticPartitionedRungeKutta.IncrementSPRK(
            computeForce: computeAccelerations,
            computeVelocity: computeVelocities,
            q0: q, p0: v, t0: t, tmax: t + duration,
            Δt: duration / Math.Ceiling(duration / timestep),
            coefficients: SymplecticPartitionedRungeKutta.Order5Optimal,
            samplingPeriod: 0, qError: qError, pError: vError, tError: tError);
      q = solution.position[0];
      v = solution.momentum[0];
      t = solution.time[0];
      updateBodies();
    }

    #region private

    private int dimension;
    private double[] q;
    private double[] qError;
    private double t;
    private double tError;
    private double[] v;
    private double[] vError;
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
        q[3 * b] = bodies[b].qx;
        q[3 * b + 1] = bodies[b].qy;
        q[3 * b + 2] = bodies[b].qz;
        v[3 * b] = bodies[b].vx;
        v[3 * b + 1] = bodies[b].vy;
        v[3 * b + 2] = bodies[b].vz;
        qError[3 * b] = bodies[b].qErrorx;
        qError[3 * b + 1] = bodies[b].qErrory;
        qError[3 * b + 2] = bodies[b].qErrorz;
        vError[3 * b] = bodies[b].vErrorx;
        vError[3 * b + 1] = bodies[b].vErrory;
        vError[3 * b + 2] = bodies[b].vErrorz;
      }
    }
    private void updateState() {
      for (int b = 0; b < bodies.Length; ++b) {
        bodies[b].qx = q[3 * b];
        bodies[b].qy = q[3 * b + 1];
        bodies[b].qz = q[3 * b + 2];
        bodies[b].vx = v[3 * b];
        bodies[b].vy = v[3 * b + 1];
        bodies[b].vz = v[3 * b + 2];
        bodies[b].qErrorx = qError[3 * b];
        bodies[b].qErrory = qError[3 * b + 1];
        bodies[b].qErrorz = qError[3 * b + 2];
        bodies[b].vErrorx = vError[3 * b];
        bodies[b].vErrory = vError[3 * b + 1];
        bodies[b].vErrorz = vError[3 * b + 2];
      }
    }

    #endregion private
  }
}