using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  public struct Scalar : IComparable, IComparable<Scalar>, IEquatable<Scalar> {
    private readonly double value;
    private Scalar(double x) { value = x; }

    #region Conversions

    public static explicit operator double(Scalar x) { return x.value; }
    public static explicit operator Scalar(double x) { return new Scalar(x); }

    #endregion Conversions

    #region Total order

    public static bool operator ==(Scalar x, Scalar y) {
      return (double)x == (double)y;
    }
    public static bool operator !=(Scalar x, Scalar y) {
      return (double)x != (double)y;
    }
    public static bool operator <(Scalar x, Scalar y) {
      return (double)x < (double)y;
    }
    public static bool operator <=(Scalar x, Scalar y) {
      return (double)x <= (double)y;
    }
    public static bool operator >(Scalar x, Scalar y) {
      return (double)x > (double)y;
    }
    public static bool operator >=(Scalar x, Scalar y) {
      return (double)x >= (double)y;
    }

    #endregion Total order

    #region Field operations

    public static Scalar operator -(Scalar x) { return (Scalar)(-(double)x); }
    public static Scalar operator +(Scalar x, Scalar y) {
      return (Scalar)((double)x + (double)y);
    }
    public static Scalar operator -(Scalar x, Scalar y) {
      return (Scalar)((double)x - (double)y);
    }
    public static Scalar operator *(Scalar x, Scalar y) {
      return (Scalar)((double)x * (double)y);
    }
    public static Scalar operator /(Scalar x, Scalar y) {
      return (Scalar)((double)x / (double)y);
    }

    #endregion Field operations

    #region Elementary functions

    public static Scalar Sin(Scalar angle) {
      return (Scalar)Math.Sin((double)angle);
    }
    public static Scalar Cos(Scalar angle) {
      return (Scalar)Math.Cos((double)angle);
    }
    public static Scalar Tan(Scalar angle) {
      return (Scalar)Math.Tan((double)angle);
    }
    public static Scalar Sqrt(Scalar angle) {
      return (Scalar)Math.Sqrt((double)angle);
    }
    public bool Equals(Scalar x) {
      return this == x;
    }
    public override int GetHashCode() {
      return this.value.GetHashCode();
    }

    #endregion Elementary functions

    #region IComparable, IComparable<Scalar> and IEquatable<Scalar> implementations

    public int CompareTo(object obj) {
      if (obj == null) {
        // MSDN: By definition, any object compares greater than (or follows)
        // null, and two null references compare equal to each other.
        return 1;
      } else if (!this.GetType().Equals(obj.GetType())) {
        throw new ArgumentException("Object is not a Scalar");
      } else {
        return this.value.CompareTo(obj);
      }
    }
    public int CompareTo(Scalar x) {
      return this.value.CompareTo(x.value);
    }
    public override bool Equals(object obj) {
      if ((obj == null) || !this.GetType().Equals(obj.GetType())) {
        return false;
      } else {
        return this == (Scalar)obj;
      }
    }

    #endregion IComparable, IComparable<Scalar> and IEquatable<Scalar> implementations
  }
}