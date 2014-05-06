using NewtonianPhysics;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Principia {
  public static class Utilities {
    // TODO(robin): Many things in this class belong somewhere else.
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
           * orbit.pos.xzy).xzy
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
           * orbit.vel.xzy).xzy
           + (orbit.referenceBody.name != "Sun"
           ? orbit.referenceBody.orbit.AbsoluteInertialVelocity()
           : new Vector3d {
             x = 10,
             y = 10,
             z = 10
           });
    }
    public static string Description<TEnum>(this TEnum source) {
      FieldInfo fieldInfo = source.GetType().GetField(source.ToString());

      DescriptionAttribute[] attributes = (DescriptionAttribute[])fieldInfo
        .GetCustomAttributes(typeof(DescriptionAttribute), false);

      if (attributes != null && attributes.Length > 0) {
        return attributes[0].Description;
      } else {
        return source.ToString();
      }
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

    public static string ToMathematica(this Vector3d v) {
      return "{" + v.x.ToString("F10") + ", "
        + v.y.ToString("F10") + ", "
        + v.z.ToString("F10") + "}";
    }
    public static string ToString(this Vector3d v,
                                  string format,
                                  bool norm = true) {
      return (norm ? "||" : "")
        + "(" + v.x.ToString(format) + ", "
        + v.y.ToString(format) + ", "
        + v.z.ToString(format) + ")"
        + (norm ? "|| = " + v.magnitude.ToString(format) : "");
    }
    public static Vector3d ToVector(this SpatialCoordinates v) {
      return new Vector3d { x = v.x, y = v.y, z = v.z };
    }
  }
}