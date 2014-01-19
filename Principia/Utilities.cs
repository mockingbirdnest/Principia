using NewtonianPhysics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Principia {
  public static class Utilities {
    public const double g0 = 9.81;
    public static SpatialCoordinates nullVector = new SpatialCoordinates {
      x = 0.0,
      y = 0.0,
      z = 0.0
    };
    // m s^-2
    public static Vector3d AbsoluteInertialPosition(this Orbit orbit) {
      double UT = Planetarium.GetUniversalTime();
      return (QuaternionD.Inverse(Planetarium.Rotation)
           * orbit.getRelativePositionAtUT(UT).xzy).xzy
        // TODO(robin): Explicit references to "Sun" EVERYWHERE.
           + (orbit.referenceBody.name != "Sun"
           ? orbit.referenceBody.orbit.AbsoluteInertialPosition()
           : new Vector3d {
             x = 10, // TODO(robin): Strange behaviour of the integrator near 0?
             y = 10,
             z = 10
           });
    }
    public static Vector3d AbsoluteInertialVelocity(this Orbit orbit) {
      double UT = Planetarium.GetUniversalTime();
      return (QuaternionD.Inverse(Planetarium.Rotation)
           * orbit.getOrbitalVelocityAtUT(UT).xzy).xzy
           + (orbit.referenceBody.name != "Sun"
           ? orbit.referenceBody.orbit.AbsoluteInertialVelocity()
           : new Vector3d {
             x = 10,
             y = 10,
             z = 10
           });
    }
    public static Body ToBody(this Vessel vessel) {
      Vector3d q = vessel.orbit.AbsoluteInertialPosition();
      Vector3d v = vessel.orbit.AbsoluteInertialVelocity();
      return new Body {
        qError = nullVector,
        q = q.ToCoordinates(),
        vError = nullVector,
        v = v.ToCoordinates(),
        gravitationalParameter = 0
      };
    }
    public static Body ToBody(this CelestialBody body) {
      Vector3d q, v;
      if (body.name == "Sun") {
        q = new Vector3d {
          x = 10,
          y = 10,
          z = 10
        };
        v = new Vector3d {
          x = 10,
          y = 10,
          z = 10
        };
      } else {
        q = body.orbit.AbsoluteInertialPosition();
        v = body.orbit.AbsoluteInertialVelocity();
      }
      return new Body {
        qError = nullVector,
        q = q.ToCoordinates(),
        vError = nullVector,
        v = v.ToCoordinates(),
        gravitationalParameter = body.GeeASL * g0 * body.Radius * body.Radius
      };
    }
    public static SpatialCoordinates ToCoordinates(this Vector3d v) {
      return new SpatialCoordinates { x = v.x, y = v.y, z = v.z };
    }
    public static Vector3d ToVector(this SpatialCoordinates v) {
      return new Vector3d { x = v.x, y = v.y, z = v.z };
    }
  }
}