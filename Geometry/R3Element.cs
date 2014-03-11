using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  public struct R3Element {
    public Scalar X, Y, Z;
    public static R3Element operator -(R3Element v) {
      return new R3Element { X = v.X, Y = v.Y, Z = v.Z };
    }
    public static R3Element operator -(R3Element left, R3Element right) {
      return new R3Element {
        X = left.X - right.X,
        Y = left.Y - right.Y,
        Z = left.Z - right.Z
      };
    }
    public static R3Element operator *(Scalar left, R3Element right) {
      return new R3Element {
        X = left * right.X,
        Y = left * right.Y,
        Z = left * right.Z
      };
    }
    public static R3Element operator *(R3Element left, Scalar right) {
      return new R3Element {
        X = left.X * right,
        Y = left.Y * right,
        Z = left.Z * right
      };
    }
    public static R3Element operator /(R3Element left, Scalar right) {
      return new R3Element {
        X = left.X / right,
        Y = left.Y / right,
        Z = left.Z / right
      };
    }
    public static R3Element operator +(R3Element left, R3Element right) {
      return new R3Element {
        X = left.X + right.X,
        Y = left.Y + right.Y,
        Z = left.Z + right.Z
      };
    }
    public R3Element Cross(R3Element right) {
      return new R3Element {
        X = this.Y * right.Z - this.Z * right.Y,
        Y = this.Z * right.X - this.X * right.Z,
        Z = this.X * right.Y - this.Y * right.X
      };
    }
    public Scalar Dot(R3Element right) {
      return this.X * right.X + this.Y * right.Y + this.Z * right.Z;
    }
  }
}