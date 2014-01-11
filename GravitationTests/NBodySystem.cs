using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace GravitationTests {

  [TestClass]
  public class NBodySystem {

    [TestMethod]
    public void CircularOrbit() {
      double[] m = { 100000, 1 };
      double[] q0 = { 0, 0, 0, 1, 0, 0 };
      double vNorm = Math.Sqrt(NewtonianGravitation.NBodySystem.G * m[0]);
      double[] v0 = { 0, 0, 0, 0, vNorm / Math.Sqrt(2), -vNorm / Math.Sqrt(2) };
      NewtonianGravitation.NBodySystem system = new NewtonianGravitation.NBodySystem(m, q0, v0);
      Integrators.SymplecticPartitionedRungeKutta.Solution simulation = system.Simulate(100000, 1, 1);
      double rmin = double.PositiveInfinity;
      double rmax = 0;
      double Δvmax = 0;
      for (int i = 0; i < simulation.Position.Length; ++i) {
        double Δq0 = simulation.Position[i][0] - simulation.Position[i][3];
        double Δq1 = simulation.Position[i][1] - simulation.Position[i][4];
        double Δq2 = simulation.Position[i][2] - simulation.Position[i][5];
        double Δv0 = simulation.Momentum[i][0] - simulation.Momentum[0][0];
        double Δv1 = simulation.Momentum[i][1] - simulation.Momentum[0][1];
        double Δv2 = simulation.Momentum[i][2] - simulation.Momentum[0][2];
        double Δv = Math.Sqrt(Δv0 * Δv0 + Δv1 * Δv1 + Δv2 * Δv2);
        double r = Math.Sqrt(Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2);
        rmin = Math.Min(r, rmin);
        rmax = Math.Max(r, rmax);
        Δvmax = Math.Max(Δv, Δvmax);
      }
      Console.WriteLine("Δvmax = " + Δvmax + ";");
      Console.WriteLine("rmin = " + rmin + ";");
      Console.WriteLine("rmax = " + rmax + ";");
      Assert.AreEqual(1, rmax);
      Assert.AreEqual(.999980000400022, rmin, 1E-15);
    }
  }
}