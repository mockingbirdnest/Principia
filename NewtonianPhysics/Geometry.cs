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
  // Structs for strong typing in a three dimensional real inner product space
  // and its Grassmann algebra, as well as an affine space.
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
  // For k-vectors v, w, the inner product (v|w) is [k]vector.InnerProduct(v, w)
  // rather than v.InnerProduct(w) so that the enclosing space is explicit.
  // Compare with the notation v.Wedge(w) for the exterior product of a k-vector
  // v and a l-vector w.
  // We choose to use operators for the vector space operations, thus keeping
  // the enclosing space implicit, i.e., v + w rather than [k]vector.Plus(v, w).
  // The left action a . v of a bivector a on v is a.ActOn(v).
  // The right action v* . a is v.ActedUponBy(a). The commutator is
  // Bivector.Commutator(a, b).

  public interface ISpace { };

  public struct Sign {
    public bool positive;
    public static explicit operator Sign(int x) {
      return new Sign { positive = x > 0 };
    }
    public static explicit operator Sign(Scalar x) {
      return new Sign { positive = x > (Scalar)0 };
    }
    public static implicit operator Scalar(Sign sgn) {
      return sgn.positive ? (Scalar)1 : (Scalar)(-1);
    }
    public static Sign operator *(Sign left, Sign right) {
      return new Sign { positive = (left.positive == right.positive) };
    }
  }

  #region Underlying weakly-typed vector

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

  #endregion Underlying weakly-typed vector

  #region Grassman algebra

  public struct Scalar : IComparable, IComparable<Scalar>, IEquatable<Scalar> {
    private double value;
    public static Scalar Cos(Scalar angle) {
      return (Scalar)Math.Cos((double)angle);
    }
    public static explicit operator double(Scalar x) { return x.value; }
    public static explicit operator Scalar(double x) {
      return new Scalar { value = x };
    }
    public static Scalar operator -(Scalar x) {
      return (Scalar)(-(double)x);
    }
    public static Scalar operator -(Scalar x, Scalar y) {
      return (Scalar)((double)x - (double)y);
    }
    public static bool operator !=(Scalar x, Scalar y) {
      return (double)x != (double)y;
    }
    public static Scalar operator *(Scalar x, Scalar y) {
      return (Scalar)((double)x * (double)y);
    }
    public static Scalar operator /(Scalar x, Scalar y) {
      return (Scalar)((double)x / (double)y);
    }
    public static Scalar operator +(Scalar x, Scalar y) {
      return (Scalar)((double)x + (double)y);
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
    public static Scalar Sin(Scalar angle) {
      return (Scalar)Math.Sin((double)angle);
    }
    public static Scalar Sqrt(Scalar angle) {
      return (Scalar)Math.Sqrt((double)angle);
    }
    public static Scalar Tan(Scalar angle) {
      return (Scalar)Math.Tan((double)angle);
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
  public struct Vector<A> where A : ISpace {
    public R3Element Coordinates;
    public static Scalar InnerProduct(Vector<A> left, Vector<A> right) {
      return left.Coordinates.Dot(right.Coordinates);
    }
    public static Vector<A> operator -(Vector<A> v) {
      return new Vector<A> { Coordinates = -v.Coordinates };
    }
    public static Vector<A> operator -(Vector<A> left, Vector<A> right) {
      return new Vector<A> {
        Coordinates = left.Coordinates - right.Coordinates
      };
    }
    public static Vector<A> operator *(Scalar left, Vector<A> right) {
      return new Vector<A> { Coordinates = left * right.Coordinates };
    }
    public static Vector<A> operator *(Vector<A> left, Scalar right) {
      return new Vector<A> { Coordinates = left.Coordinates * right };
    }
    public static Vector<A> operator /(Vector<A> left, Scalar right) {
      return new Vector<A> { Coordinates = left.Coordinates / right };
    }
    public static Vector<A> operator +(Vector<A> left, Vector<A> right) {
      return new Vector<A> {
        Coordinates = left.Coordinates + right.Coordinates
      };
    }
    public static Vector<A> ToFrom(Point<A> left, Point<A> right) {
      return new Vector<A> {
        Coordinates = left.Coordinates - right.Coordinates
      };
    }
    public Vector<A> ActedUponBy(BiVector<A> right) {
      return new Vector<A> {
        Coordinates = this.Coordinates.Cross(right.Coordinates)
      };
    }
    public Point<A> Translate(Point<A> right) {
      return new Point<A> {
        Coordinates = this.Coordinates + right.Coordinates
      };
    }
    public BiVector<A> Wedge(Vector<A> right) {
      return new BiVector<A> {
        Coordinates = this.Coordinates.Cross(right.Coordinates)
      };
    }
    public TriVector<A> Wedge(BiVector<A> right) {
      return new TriVector<A> {
        Coordinate = this.Coordinates.Dot(right.Coordinates)
      };
    }
  }
  public struct BiVector<A> where A : ISpace {
    public R3Element Coordinates;
    public static BiVector<A> Commutator(BiVector<A> left, BiVector<A> right) {
      return new BiVector<A> {
        Coordinates = left.Coordinates.Cross(right.Coordinates)
      };
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
      return new BiVector<A> { Coordinates = -v.Coordinates };
    }
    public static BiVector<A> operator -(BiVector<A> left, BiVector<A> right) {
      return new BiVector<A> {
        Coordinates = left.Coordinates - right.Coordinates
      };
    }
    public static BiVector<A> operator *(Scalar left, BiVector<A> right) {
      return new BiVector<A> { Coordinates = left * right.Coordinates };
    }
    public static BiVector<A> operator *(BiVector<A> left, Scalar right) {
      return new BiVector<A> { Coordinates = left.Coordinates * right };
    }
    public static BiVector<A> operator /(BiVector<A> left, Scalar right) {
      return new BiVector<A> { Coordinates = left.Coordinates / right };
    }
    public static BiVector<A> operator +(BiVector<A> left, BiVector<A> right) {
      return new BiVector<A> {
        Coordinates = left.Coordinates + right.Coordinates
      };
    }
    public Vector<A> ActOn(Vector<A> right) {
      return new Vector<A> {
        Coordinates = this.Coordinates.Cross(right.Coordinates)
      };
    }
    public TriVector<A> Wedge(Vector<A> right) {
      return new TriVector<A> {
        Coordinate = this.Coordinates.Dot(right.Coordinates)
      };
    }
  }
  public struct TriVector<A> where A : ISpace {
    public Scalar Coordinate;
  }

  #endregion Grassman algebra

  #region Affine space

  public struct Point<A> where A : ISpace {
    public R3Element Coordinates;
    // A convex combination of positions is a position.
    public static Point<A> Barycenter(Point<A> q1, Scalar λ1,
                                      Point<A> q2, Scalar λ2) {
      return new Point<A> {
        Coordinates = (q1.Coordinates * λ1 + q2.Coordinates * λ2) / (λ1 + λ2)
      };
    }
  }

  #endregion Affine space

  #region Maps between Grassmann algebras

  // A permutation of the coordinates. Obviously not coordinate-free, but
  // practical. There are no precision losses when composing or applying a
  // coordinate permutation.
  public struct Permutation<A, B>
    where A : ISpace
    where B : ISpace {
    public CoordinatePermutation BasisImage;
    public enum CoordinatePermutation {
      XYZ = 1, YZX = 2, ZXY = 3,
      XZY = -1, ZYX = -2, YXZ = -3
    }
    public Sign Determinant { get { return (Sign)(int)BasisImage; } }
    public static R3Element operator *(Permutation<A, B> left,
                                       R3Element right) {
      switch (left.BasisImage) {
        case CoordinatePermutation.XYZ:
          return right;
        case CoordinatePermutation.XZY:
          return new R3Element {
            X = right.X,
            Y = right.Z,
            Z = right.Y
          };
        case CoordinatePermutation.YXZ:
          return new R3Element {
            X = right.Y,
            Y = right.X,
            Z = right.Z
          };
        case CoordinatePermutation.YZX:
          return new R3Element {
            X = right.Y,
            Y = right.Z,
            Z = right.X
          };
        case CoordinatePermutation.ZXY:
          return new R3Element {
            X = right.Z,
            Y = right.X,
            Z = right.Y
          };
        case CoordinatePermutation.ZYX:
          return new R3Element {
            X = right.Z,
            Y = right.Y,
            Z = right.X
          };
        default:
          // Stupid language.
          Console.WriteLine("CoordinatePermutation.BasisImage "
                            + "was out of scope.");
          return new R3Element { };
      }
    }
    public OrthogonalTransformation<A, B> Forget() {
      Rotation<A, B> specialOrthogonalMap;
      switch (BasisImage) {
        case CoordinatePermutation.XYZ:
          specialOrthogonalMap.RealPart = (Scalar)1;
          specialOrthogonalMap.ImaginaryPart = new R3Element {
            X = (Scalar)0,
            Y = (Scalar)0,
            Z = (Scalar)0
          };
          break;
        case CoordinatePermutation.YZX:
          specialOrthogonalMap.RealPart = (Scalar)(-.5);
          specialOrthogonalMap.ImaginaryPart = new R3Element {
            X = (Scalar).5,
            Y = (Scalar).5,
            Z = (Scalar).5
          };
          break;
        case CoordinatePermutation.ZXY:
          specialOrthogonalMap.RealPart = (Scalar).5;
          specialOrthogonalMap.ImaginaryPart = new R3Element {
            X = (Scalar).5,
            Y = (Scalar).5,
            Z = (Scalar).5
          };
          break;
        case CoordinatePermutation.XZY:
          specialOrthogonalMap.RealPart = (Scalar)0;
          specialOrthogonalMap.ImaginaryPart = new R3Element {
            X = (Scalar)0,
            Y = -Scalar.Sqrt((Scalar)2) / (Scalar)2,
            Z = Scalar.Sqrt((Scalar)2) / (Scalar)2
          };
          break;
        case CoordinatePermutation.YXZ:
          specialOrthogonalMap.RealPart = (Scalar)0;
          specialOrthogonalMap.ImaginaryPart = new R3Element {
            X = -Scalar.Sqrt((Scalar)2) / (Scalar)2,
            Y = Scalar.Sqrt((Scalar)2) / (Scalar)2,
            Z = (Scalar)0
          };
          break;
        case CoordinatePermutation.ZYX:
          specialOrthogonalMap.RealPart = (Scalar)0;
          specialOrthogonalMap.ImaginaryPart = new R3Element {
            X = -Scalar.Sqrt((Scalar)2) / (Scalar)2,
            Y = (Scalar)0,
            Z = Scalar.Sqrt((Scalar)2) / (Scalar)2
          };
          break;
        default:
          // Stupid language.
          Console.WriteLine("CoordinatePermutation.BasisImage "
                            + "was out of scope.");
          return new OrthogonalTransformation<A, B>();
      }
      return new OrthogonalTransformation<A, B> {
        Determinant = this.Determinant,
        SpecialOrthogonalMap = specialOrthogonalMap
      };
    }
    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B> { Coordinates = this * right.Coordinates };
    }
    public BiVector<A> ActOn(BiVector<A> right) {
      return new BiVector<A> {
        Coordinates = this.Determinant * (this * right.Coordinates)
      };
    }
    public TriVector<A> ActOn(TriVector<A> right) {
      return new TriVector<A> {
        Coordinate = this.Determinant * right.Coordinate
      };
    }
  }

  // The rotation is modeled as a quaternion.
  public struct Rotation<A, B>
    where A : ISpace
    where B : ISpace {
    public Scalar RealPart;
    public R3Element ImaginaryPart;
    public Rotation<B, A> Inverse {
      get {
        return new Rotation<B, A> {
          RealPart = RealPart,
          ImaginaryPart = -ImaginaryPart
        };
      }
    }
    public static R3Element operator *(Rotation<A, B> left,
                                       R3Element right) {
      // TODO(egg): Optimise.
      return Maps.Compose<B, A, B>(
        Maps.Compose<A, A, B>(
          left,
          new Rotation<A, A> { RealPart = (Scalar)0, ImaginaryPart = right }),
        left.Inverse).ImaginaryPart;
    }
    public OrthogonalTransformation<A, B> Forget() {
      return new OrthogonalTransformation<A, B> {
        Determinant = (Sign)1,
        SpecialOrthogonalMap = this
      };
    }
    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B> { Coordinates = this * right.Coordinates };
    }
  }
  public struct OrthogonalTransformation<A, B>
    where A : ISpace
    where B : ISpace {
    // The orthogonal transformation is modeled as a rotoinversion.
    public Sign Determinant;
    public Rotation<A, B> SpecialOrthogonalMap;

    public static R3Element operator *(OrthogonalTransformation<A, B> left,
                                      R3Element right) {
      return left.SpecialOrthogonalMap * (left.Determinant * right);
    }
    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B> { Coordinates = this * right.Coordinates };
    }
    public BiVector<A> ActOn(BiVector<A> right) {
      return new BiVector<A> {
        Coordinates = this.SpecialOrthogonalMap * right.Coordinates
      };
    }
    public TriVector<A> ActOn(TriVector<A> right) {
      return new TriVector<A> {
        Coordinate = this.Determinant * right.Coordinate
      };
    }
  }

  #endregion Maps between Grassmann algebras

  #region Maps between affine spaces

  public struct RigidTransformation<A, B>
    where A : ISpace
    where B : ISpace {
    public Rotation<A, B> orthogonalMap;
    public Vector<B> translation;
  }
  public struct EuclideanTransformation<A, B>
    where A : ISpace
    where B : ISpace {
    public OrthogonalTransformation<A, B> orthogonalMap;
    public Vector<B> translation;
  }

  #endregion Maps between affine spaces

  #region Composition of maps

  public static class Maps {
    public static Rotation<A, C> Compose<A, B, C>(Rotation<B, C> left,
                                                  Rotation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new Rotation<A, C> {
        RealPart = left.RealPart * right.RealPart
                   - left.ImaginaryPart.Dot(right.ImaginaryPart),
        ImaginaryPart = left.RealPart * right.ImaginaryPart
                        + right.RealPart * left.ImaginaryPart
                        + left.ImaginaryPart.Cross(right.ImaginaryPart)
      };
    }
    public static OrthogonalTransformation<A, C> Compose<A, B, C>(
      OrthogonalTransformation<B, C> left,
      OrthogonalTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new OrthogonalTransformation<A, C> {
        Determinant = left.Determinant * right.Determinant,
        SpecialOrthogonalMap = Compose<A, B, C>(left.SpecialOrthogonalMap,
                                                right.SpecialOrthogonalMap)
      };
    }
    public static EuclideanTransformation<A, C> Compose<A, B, C>(
      EuclideanTransformation<B, C> left,
      EuclideanTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new EuclideanTransformation<A, C> {
        orthogonalMap = Compose<A, B, C>(left.orthogonalMap,
                                         right.orthogonalMap),
        translation = left.translation
                      + left.orthogonalMap.ActOn(right.translation)
      };
    }
    public static RigidTransformation<A, C> Compose<A, B, C>(
      RigidTransformation<B, C> left,
      RigidTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new RigidTransformation<A, C> {
        orthogonalMap = Compose<A, B, C>(left.orthogonalMap,
                                         right.orthogonalMap),
        translation = left.translation
                      + left.orthogonalMap.ActOn(right.translation)
      };
    }
  }

  #endregion Composition of maps
}