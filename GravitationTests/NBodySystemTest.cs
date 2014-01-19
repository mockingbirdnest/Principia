using Microsoft.VisualStudio.TestTools.UnitTesting;
using NewtonianPhysics;
using System;

namespace GravitationTests {
  [TestClass]
  public class NBodySystemTest {
    [TestMethod]
    public void CircularOrbit() {
      double m0 = 100000;
      double m1 = 1;
      double vNorm = Math.Sqrt(NBodySystem.G * m0);
      Body[] bodies = new Body[] {
        new Body { gravitationalParameter = m0 * NBodySystem.G },
        new Body { q= new SpatialCoordinates{x = 1},
          gravitationalParameter = m1 * NBodySystem.G ,
          v=new SpatialCoordinates{y= vNorm / Math.Sqrt(2) ,z = -vNorm / Math.Sqrt(2)} }};
      double[] q0 = { 0, 0, 0, 1, 0, 0 };
      // r = 1, v0 = 0, v1^2 =  G m0, so U = - G m0 m1, T = 1/2 G m0 m1
      double Estart = -(NBodySystem.G * m0 * m1) / 2;
      NBodySystem system = new NBodySystem(bodies, 0);
      system.Evolve(100000, 1);
      double Δqx = bodies[0].q.x - bodies[1].q.x;
      double Δqy = bodies[0].q.y - bodies[1].q.y;
      double Δqz = bodies[0].q.z - bodies[1].q.z;
      double r = Math.Sqrt(Δqx * Δqx + Δqy * +Δqy + Δqz * Δqz);
      double v0Squared = bodies[0].v.x * bodies[0].v.x
                       + bodies[0].v.y * bodies[0].v.y
                       + bodies[0].v.z * bodies[0].v.z;
      double v1Squared = bodies[1].v.x * bodies[1].v.x
                       + bodies[1].v.y * bodies[1].v.y
                       + bodies[1].v.z * bodies[1].v.z;
      double Eend = -NBodySystem.G * m0 * m1 / r
                  + (m0 * v0Squared + m1 * v1Squared) / 2;
      double Eerror = Math.Abs((Eend - Estart) / Estart);
      Console.WriteLine("v0^2 = " + v0Squared + ";");
      Console.WriteLine("v1^2 = " + v1Squared + ";");
      Console.WriteLine("r = " + r + ";");
      Console.WriteLine("ΔE/E0 = " + Eerror + ";");
      Assert.IsTrue(Eerror < 1E-15);
      Assert.IsTrue(Eerror > 1E-16);
    }
  }
}