using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NewtonianPhysics {
  public struct Body {
    // We use the gravitational parameter μ = G M in order not to accumulate
    // unit roundoffs from repeated multiplications by G. Note that in KSP, the
    // gravitational parameter is computed from the mass as G M, but the mass is
    // itself computed from the radius and acceleration due to gravity at sea
    // level as M = g0 r^2 / G. This is silly (and introduces an---admittedly
    // tiny---error), so the gravitational parameter should ideally be computed
    // by the user as μ = g0 r^2.
    public double gravitationalParameter;

    // Errors from compensated summation.
    public double qErrorx, qErrory, qErrorz, vErrorx, vErrory, vErrorz;

    // We don't use KSP's Vector3d for the following reasons:
    // 1. Sloppy numerics (there are places which explicitly underflow values
    // below 1E-6, angles are constantly converted to degrees, etc.)
    // 2. We want the NewtonianPhysics assembly to be independent from Unity and
    // KSP (the Simulator is a standalone application and should not require KSP
    // to run).
    // We don't use double arrays because they are reference types (silly
    // language) and the number of dimensions is not likely to change anyway.
    public double qx, qy, qz, vx, vy, vz;
    public readonly bool massless {
      get { return gravitationalParameter == 0; }
    }
  }
}