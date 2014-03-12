using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
// Strong typing between Grassmann algebras (e.g., between reference frames) is
// enforced by genericity over tags deriving from the empty interface ISpace.
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
// For k-vectors v, w over the underlying space V, the inner product (v|w) is
// [k]vector<V>.InnerProduct(v, w) rather than v.InnerProduct(w) so that the
// enclosing space is explicit. Compare with the notation v.Wedge(w) for the
// exterior product of a k-vector v and a l-vector w.
// We choose to use operators for the vector space operations, thus keeping
// the enclosing space implicit, i.e., v + w rather than
// [k]vector<V>.Plus(v, w). The left action a . v of a bivector a on v is
// a.ActOn(v). The right action v* . a is v.ActedUponBy(a). The commutator is
// Bivector<V>.Commutator(a, b).

namespace Geometry {
  public struct Vector<A> where A : ISpace {
    public readonly R3Element Coordinates;
    public Vector(R3Element coordinates) {
      Coordinates = coordinates;
    }
    public static Scalar InnerProduct(Vector<A> left, Vector<A> right) {
      return left.Coordinates.Dot(right.Coordinates);
    }
    public static Vector<A> operator -(Vector<A> v) {
      return new Vector<A>(-v.Coordinates);
    }
    public static Vector<A> operator -(Vector<A> left, Vector<A> right) {
      return new Vector<A>(left.Coordinates - right.Coordinates);
    }
    public static Vector<A> operator *(Scalar left, Vector<A> right) {
      return new Vector<A>(left * right.Coordinates);
    }
    public static Vector<A> operator *(Vector<A> left, Scalar right) {
      return new Vector<A>(left.Coordinates * right);
    }
    public static Vector<A> operator /(Vector<A> left, Scalar right) {
      return new Vector<A>(left.Coordinates / right);
    }
    public static Vector<A> operator +(Vector<A> left, Vector<A> right) {
      return new Vector<A>(left.Coordinates + right.Coordinates);
    }
    public static Vector<A> ToFrom(Point<A> to, Point<A> from) {
      return new Vector<A>(to.Coordinates - from.Coordinates);
    }
    public Vector<A> ActedUponBy(BiVector<A> right) {
      return new Vector<A>(this.Coordinates.Cross(right.Coordinates));
    }
    public Point<A> Translate(Point<A> right) {
      return new Point<A>(this.Coordinates + right.Coordinates);
    }
    public BiVector<A> Wedge(Vector<A> right) {
      return new BiVector<A>(this.Coordinates.Cross(right.Coordinates));
    }
    public TriVector<A> Wedge(BiVector<A> right) {
      return new TriVector<A>(this.Coordinates.Dot(right.Coordinates));
    }
  }
  public struct BiVector<A> where A : ISpace {
    public readonly R3Element Coordinates;
    public BiVector(R3Element coordinates) {
      Coordinates = coordinates;
    }
    public static BiVector<A> Commutator(BiVector<A> left, BiVector<A> right) {
      return new BiVector<A>(left.Coordinates.Cross(right.Coordinates));
    }
    public static Rotation<A, A> Exp(BiVector<A> infinitesimalRotation) {
      Scalar angle = Scalar.Sqrt(BiVector<A>.InnerProduct(infinitesimalRotation,
                                                       infinitesimalRotation));
      return new Rotation<A, A> {
        RealPart = Scalar.Cos(angle / (Scalar)2),
        ImaginaryPart = infinitesimalRotation.Coordinates / angle
                        * Scalar.Sin(angle / (Scalar)2)
      };
    }
    public static Scalar InnerProduct(BiVector<A> left, BiVector<A> right) {
      return left.Coordinates.Dot(right.Coordinates);
    }
    public static BiVector<A> operator -(BiVector<A> v) {
      return new BiVector<A>(-v.Coordinates);
    }
    public static BiVector<A> operator -(BiVector<A> left, BiVector<A> right) {
      return new BiVector<A>(left.Coordinates - right.Coordinates);
    }
    public static BiVector<A> operator *(Scalar left, BiVector<A> right) {
      return new BiVector<A>(left * right.Coordinates);
    }
    public static BiVector<A> operator *(BiVector<A> left, Scalar right) {
      return new BiVector<A>(left.Coordinates * right);
    }
    public static BiVector<A> operator /(BiVector<A> left, Scalar right) {
      return new BiVector<A>(left.Coordinates / right);
    }
    public static BiVector<A> operator +(BiVector<A> left, BiVector<A> right) {
      return new BiVector<A>(left.Coordinates + right.Coordinates);
    }
    public Vector<A> ActOn(Vector<A> right) {
      return new Vector<A>(this.Coordinates.Cross(right.Coordinates));
    }
    public TriVector<A> Wedge(Vector<A> right) {
      return new TriVector<A>(this.Coordinates.Dot(right.Coordinates));
    }
  }
  public struct TriVector<A> where A : ISpace {
    public readonly Scalar Coordinate;
    public TriVector(Scalar coordinate) {
      Coordinate = coordinate;
    }
  }
}