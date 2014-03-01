using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

// This namespace is a monstrosity as far as code replication is concerned, but
// there is no way around that short of rewriting the whole thing in C++/CLI.
// This might be done eventually (especially if we need to interface with native
// C++), but is not a high priority. It would make the management of physical
// quantities easier and more systematic.
namespace NewtonianPhysics.Geometry {

  #region Linear Algebra

  // Structs for strong typing in a three dimensional real inner product space
  // and its Grassmann algebra.
  // Double-precision floating point numbers, representing reals, are wrapped in
  // the struct Scalar to ensure strong typing. Three-dimensional data is stored
  // as R3Element.
  // We use the Grassman algebra in order to avoid the confusion between pseudo-
  // vectors and vectors and between pseudoscalars and scalars. This allows for
  // handling of indirect changes of coordinates. As a result, there is no cross
  // product except for the underlying data type R3Element.
  // We treat the isomorphism between a space and its dual as implicit, as we
  // are dealing with inner product spaces.
  // TODO(robin): We may want to change this when we get to rotating-pulsating
  // reference frames.
  // A reminder of the equivalents in the Grassman algebra of the cross
  // product; We identify V ^ V with the Lie algebra so(V).
  // . pseudovector = vector x vector, as in L = r x p, becomes the wedge
  //   product,
  //   bivector = vector ^ vector, L = r ^ p.
  // . vector = pseudovector x vector, as in a_Coriolis = 2 Ω x v, becomes
  //   the left action of the bivector,
  //   vector = bivector . vector,  a_Coriolis = 2 Ω . v.
  // . vector = vector x pseudovector, as in F = q v x B, becomes
  //   the dual action,
  //   vector = bivector . vector,  F* = q v* . B.
  // . pseudovector = pseudovector x pseudovector, as in the composition of
  //   rotation vectors Θ12 = α Θ1 + β Θ2 + γ Θ1 x Θ2, becomes the Lie bracket
  //   of so(V).
  //   bivector = [bivector, bivector], Θ12 = α Θ1 + β Θ2 + γ [Θ1, Θ2].
  // . pseudoscalar = (pseudovector|vector), as in Φ = (B|S), becomes the wedge
  //   product,
  //   trivector = bivector ^ vector, Φ = B ^ S.
  // NOTATION:
  // For k-vectors v, w, the inner product (v|w) is [k]Vector.InnerProduct(v, w)
  // rather than v.InnerProduct(w) so that the enclosing space is explicit.
  // Compare with the notation v.Wedge(w) for the exterior product of a k-vector
  // v and a l-vector w.
  // We choose to use operators for the vector space operations, thus keeping
  // the enclosing space implicit, i.e., v + w rather than [k]Vector.Plus(v, w).
  // The left action a . v of a bivector a on v is a.ActOn(v).
  // The right action v* . a is v.ActedUpon(a). The commutator is
  // Bivector.Commutator(a, b).

  public struct Bivector {
    public R3Element coordinates;
    public static Bivector Commutator(Bivector left, Bivector right) {
      return new Bivector {
        coordinates = left.coordinates.Cross(right.coordinates)
      };
    }
    public static Scalar InnerProduct(Bivector left, Bivector right) {
      return left.coordinates.Dot(right.coordinates);
    }
    public static Bivector operator -(Bivector v) {
      return new Bivector { coordinates = -v.coordinates };
    }
    public static Bivector operator -(Bivector left, Bivector right) {
      return new Bivector { coordinates = left.coordinates - right.coordinates };
    }
    public static Bivector operator *(Scalar left, Bivector right) {
      return new Bivector { coordinates = left * right.coordinates };
    }
    public static Bivector operator *(Bivector left, Scalar right) {
      return new Bivector { coordinates = left.coordinates * right };
    }
    public static Bivector operator /(Bivector left, Scalar right) {
      return new Bivector { coordinates = left.coordinates / right };
    }
    public static Bivector operator +(Bivector left, Bivector right) {
      return new Bivector { coordinates = left.coordinates + right.coordinates };
    }
    public Vector ActOn(Vector right) {
      return new Vector {
        coordinates = this.coordinates.Cross(right.coordinates)
      };
    }
    public TriVector Wedge(Vector right) {
      return new TriVector {
        coordinate = this.coordinates.Dot(right.coordinates)
      };
    }
  }
  public struct R3Element {
    public Scalar x, y, z;
    public static R3Element operator -(R3Element v) {
      return new R3Element { x = v.x, y = v.y, z = v.z };
    }
    public static R3Element operator -(R3Element left, R3Element right) {
      return new R3Element {
        x = left.x - right.x,
        y = left.y - right.y,
        z = left.z - right.z
      };
    }
    public static R3Element operator *(Scalar left, R3Element right) {
      return new R3Element {
        x = left * right.x,
        y = left * right.y,
        z = left * right.z
      };
    }
    public static R3Element operator *(R3Element left, Scalar right) {
      return new R3Element {
        x = left.x * right,
        y = left.y * right,
        z = left.z * right
      };
    }
    public static R3Element operator /(R3Element left, Scalar right) {
      return new R3Element {
        x = left.x / right,
        y = left.y / right,
        z = left.z / right
      };
    }
    public static R3Element operator +(R3Element left, R3Element right) {
      return new R3Element {
        x = left.x + right.x,
        y = left.y + right.y,
        z = left.z + right.z
      };
    }
    public R3Element Cross(R3Element right) {
      return new R3Element {
        x = this.y * right.z - this.z * right.y,
        y = this.z * right.x - this.x * right.z,
        z = this.x * right.y - this.y * right.x
      };
    }
    public Scalar Dot(R3Element right) {
      return this.x * right.x + this.y * right.y + this.z * right.z;
    }
  }
  public struct Scalar : IComparable, IComparable<Scalar>, IEquatable<Scalar> {
    private double value;
    public static explicit operator double(Scalar x) { return x.value; }
    public static implicit operator Scalar(double x) {
      return new Scalar { value = x };
    }
    public static Scalar operator -(Scalar x, Scalar y) {
      return (double)x - (double)y;
    }
    public static bool operator !=(Scalar x, Scalar y) {
      return (double)x != (double)y;
    }
    public static Scalar operator *(Scalar x, Scalar y) {
      return (double)x * (double)y;
    }
    public static Scalar operator /(Scalar x, Scalar y) {
      return (double)x / (double)y;
    }
    public static Scalar operator +(Scalar x, Scalar y) {
      return (double)x + (double)y;
    }
    public static bool operator <(Scalar x, Scalar y) {
      return (double)x < (double)y;
    }
    public static bool operator <=(Scalar x, Scalar y) {
      return (double)x <= (double)y;
    }
    public static bool operator ==(Scalar x, Scalar y) {
      return (double)x == (double)y;
    }
    public static bool operator >(Scalar x, Scalar y) {
      return (double)x > (double)y;
    }
    public static bool operator >=(Scalar x, Scalar y) {
      return (double)x >= (double)y;
    }
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
    public bool Equals(Scalar x) {
      return this == x;
    }
    public override int GetHashCode() {
      return this.value.GetHashCode();
    }
  }
  public struct TriVector {
    public Scalar coordinate;
  }
  public struct Vector {
    public R3Element coordinates;
    public static Scalar InnerProduct(Vector left, Vector right) {
      return left.coordinates.Dot(right.coordinates);
    }
    public static Vector operator -(Vector v) {
      return new Vector { coordinates = -v.coordinates };
    }
    public static Vector operator -(Vector left, Vector right) {
      return new Vector { coordinates = left.coordinates - right.coordinates };
    }
    public static Vector operator *(Scalar left, Vector right) {
      return new Vector { coordinates = left * right.coordinates };
    }
    public static Vector operator *(Vector left, Scalar right) {
      return new Vector { coordinates = left.coordinates * right };
    }
    public static Vector operator /(Vector left, Scalar right) {
      return new Vector { coordinates = left.coordinates / right };
    }
    public static Vector operator +(Vector left, Vector right) {
      return new Vector { coordinates = left.coordinates + right.coordinates };
    }
    public Vector ActedUpon(Bivector right) {
      return new Vector {
        coordinates = this.coordinates.Cross(right.coordinates)
      };
    }
    public Bivector Wedge(Vector right) {
      return new Bivector {
        coordinates = this.coordinates.Cross(right.coordinates)
      };
    }
    public TriVector Wedge(Bivector right) {
      return new TriVector {
        coordinate = this.coordinates.Dot(right.coordinates)
      };
    }
  }

  #endregion Linear Algebra

  public interface IReferenceFrame { };
}