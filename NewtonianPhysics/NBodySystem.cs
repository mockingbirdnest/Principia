using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NewtonianPhysics {

  public class NBodySystem {

    public NBodySystem(double[] masses, double[] q0, double[] v0) {
      this.masses = masses;
      this.v0 = v0;
      this.q0 = q0;
      this.bodies = masses.Length;
      this.dimension = 3 * bodies;
    }

    public Integrators.SymplecticPartitionedRungeKutta.Solution Simulate(double tmax, double Δt, int samplingPeriod) {
      return Integrators.SymplecticPartitionedRungeKutta.IncrementSPRK(
        computeAccelerations, computeVelocities,
        q0, v0,
        0, tmax, Δt,
        Integrators.SymplecticPartitionedRungeKutta.Order5Coefficients,
        samplingPeriod);
    }

    private void computeAccelerations(double[] q, double t, ref double[] result) {
      for (int k = 0; k < dimension; ++k) {
        result[k] = 0;
      }
      for (int b1 = 0; b1 < bodies; ++b1) {
        for (int b2 = b1 + 1; b2 < bodies; ++b2) {
          double Δq0 = q[3 * b1] - q[3 * b2];
          double Δq1 = q[3 * b1 + 1] - q[3 * b2 + 1];
          double Δq2 = q[3 * b1 + 2] - q[3 * b2 + 2];
          double squaredDistance = Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
          double denominator = squaredDistance * Math.Sqrt(squaredDistance);
          double m1OverRSquared = masses[b1] / denominator;
          double m2OverRSquared = masses[b2] / denominator;
          result[3 * b1] -= Δq0 * m2OverRSquared;
          result[3 * b1 + 1] -= Δq1 * m2OverRSquared;
          result[3 * b1 + 2] -= Δq2 * m2OverRSquared;
          // Lex. III. Actioni contrariam semper & æqualem esse reactionem: sive
          // corporum duorum actiones in se mutuo semper esse æquales & in
          // partes contrarias dirigi.
          result[3 * b2] += Δq0 * m1OverRSquared;
          result[3 * b2 + 1] += Δq1 * m1OverRSquared;
          result[3 * b2 + 2] += Δq2 * m1OverRSquared;
        }
      }

      for (int k = 0; k < dimension; ++k) {
        result[k] *= G;
      }
    }

    private void computeVelocities(double[] v, ref double[] result) {
      result = v;
    }

    public readonly double[] masses;

    public readonly double[] q0;

    public readonly double[] v0;

    public readonly int bodies;
    public readonly int dimension;
    public const double G = 6.67384E-11;
  }
}