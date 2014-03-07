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
  // Vector<V>s and Vector<V>s and between pseudoscalars and scalars. This allows for
  // handling of indirect changes of coordinates. As a result, there is no cross
  // product except for the underlying data type R3Element.
  // We treat the isomorphism between a space and its dual as implicit, as we
  // are dealing with inner product spaces.
  // TODO(robin): We may want to change this when we get to rotating-pulsating
  // reference frames.
  // A reminder of the equivalents in the Grassman algebra of the cross
  // product; We identify V ^ V with the Lie algebra so(V).
  // . pseudoVector<V> = Vector<V> x Vector<V>, as in L = r x p, becomes the wedge
  //   product,
  //   biVector<V> = Vector<V> ^ Vector<V>, L = r ^ p.
  // . Vector<V> = pseudoVector<V> x Vector<V>, as in a_Coriolis = 2 Ω x v, becomes
  //   the left action of the biVector<V>,
  //   Vector<V> = biVector<V> . Vector<V>,  a_Coriolis = 2 Ω . v.
  // . Vector<V> = Vector<V> x pseudoVector<V>, as in F = q v x B, becomes
  //   the dual action,
  //   Vector<V> = biVector<V> . Vector<V>,  F* = q v* . B.
  // . pseudoVector<V> = pseudoVector<V> x pseudoVector<V>, as in the composition of
  //   rotation Vector<V>s Θ12 = α Θ1 + β Θ2 + γ Θ1 x Θ2, becomes the Lie bracket
  //   of so(V).
  //   biVector<V> = [biVector<V>, biVector<V>], Θ12 = α Θ1 + β Θ2 + γ [Θ1, Θ2].
  // . pseudoscalar = (pseudoVector<V>|Vector<V>), as in Φ = (B|S), becomes the wedge
  //   product,
  //   triVector<V> = biVector<V> ^ Vector<V>, Φ = B ^ S.
  // NOTATION:
  // For k-Vector<V>s v, w, the inner product (v|w) is [k]Vector<V>.InnerProduct(v, w)
  // rather than v.InnerProduct(w) so that the enclosing space is explicit.
  // Compare with the notation v.Wedge(w) for the exterior product of a k-Vector<V>
  // v and a l-Vector<V> w.
  // We choose to use operators for the Vector<V> space operations, thus keeping
  // the enclosing space implicit, i.e., v + w rather than [k]Vector<V>.Plus(v, w).
  // The left action a . v of a biVector<V> a on v is a.ActOn(v).
  // The right action v* . a is v.ActedUponBy(a). The commutator is
  // BiVector<V>.Commutator(a, b).

  public interface ISpace { };

  public struct Sign {
    public bool positive;
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

  #region Underlying weakly-typed Vector<V>

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

  #endregion Underlying weakly-typed Vector<V>

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
  public struct Vector<V> where V : ISpace {
    public R3Element coordinates;
    public static Scalar InnerProduct(Vector<V> left, Vector<V> right) {
      return left.coordinates.Dot(right.coordinates);
    }
    public static Vector<V> operator -(Vector<V> v) {
      return new Vector<V> { coordinates = -v.coordinates };
    }
    public static Vector<V> operator -(Vector<V> left, Vector<V> right) {
      return new Vector<V> { coordinates = left.coordinates - right.coordinates };
    }
    public static Vector<V> operator *(Scalar left, Vector<V> right) {
      return new Vector<V> { coordinates = left * right.coordinates };
    }
    public static Vector<V> operator *(Vector<V> left, Scalar right) {
      return new Vector<V> { coordinates = left.coordinates * right };
    }
    public static Vector<V> operator /(Vector<V> left, Scalar right) {
      return new Vector<V> { coordinates = left.coordinates / right };
    }
    public static Vector<V> operator +(Vector<V> left, Vector<V> right) {
      return new Vector<V> { coordinates = left.coordinates + right.coordinates };
    }
    public static Vector<V> ToFrom(Point<V> left, Point<V> right) {
      return new Vector<V> { coordinates = left.coordinates - right.coordinates };
    }
    public Vector<V> ActedUponBy(BiVector<V> right) {
      return new Vector<V> {
        coordinates = this.coordinates.Cross(right.coordinates)
      };
    }
    public Point<V> Translate(Point<V> right) {
      return new Point<V> {
        coordinates = this.coordinates + right.coordinates
      };
    }
    public BiVector<V> Wedge(Vector<V> right) {
      return new BiVector<V> {
        coordinates = this.coordinates.Cross(right.coordinates)
      };
    }
    public TriVector<V> Wedge(BiVector<V> right) {
      return new TriVector<V> {
        coordinate = this.coordinates.Dot(right.coordinates)
      };
    }
  }
  public struct BiVector<V> where V : ISpace {
    public R3Element coordinates;
    public static BiVector<V> Commutator(BiVector<V> left, BiVector<V> right) {
      return new BiVector<V> {
        coordinates = left.coordinates.Cross(right.coordinates)
      };
    }
    public static Rotation<V, V> Exp(BiVector<V> infinitesimalRotation) {
      Scalar angle = Scalar.Sqrt(BiVector<V>.InnerProduct(infinitesimalRotation,
                                                       infinitesimalRotation));
      return new Rotation<V, V> {
        realPart = Scalar.Cos(angle / (Scalar)2),
        imaginaryPart = infinitesimalRotation.coordinates / angle
                        * Scalar.Sin(angle / (Scalar)2)
      };
    }
    public static Scalar InnerProduct(BiVector<V> left, BiVector<V> right) {
      return left.coordinates.Dot(right.coordinates);
    }
    public static BiVector<V> operator -(BiVector<V> v) {
      return new BiVector<V> { coordinates = -v.coordinates };
    }
    public static BiVector<V> operator -(BiVector<V> left, BiVector<V> right) {
      return new BiVector<V> {
        coordinates = left.coordinates - right.coordinates
      };
    }
    public static BiVector<V> operator *(Scalar left, BiVector<V> right) {
      return new BiVector<V> { coordinates = left * right.coordinates };
    }
    public static BiVector<V> operator *(BiVector<V> left, Scalar right) {
      return new BiVector<V> { coordinates = left.coordinates * right };
    }
    public static BiVector<V> operator /(BiVector<V> left, Scalar right) {
      return new BiVector<V> { coordinates = left.coordinates / right };
    }
    public static BiVector<V> operator +(BiVector<V> left, BiVector<V> right) {
      return new BiVector<V> {
        coordinates = left.coordinates + right.coordinates
      };
    }
    public Vector<V> ActOn(Vector<V> right) {
      return new Vector<V> {
        coordinates = this.coordinates.Cross(right.coordinates)
      };
    }
    public TriVector<V> Wedge(Vector<V> right) {
      return new TriVector<V> {
        coordinate = this.coordinates.Dot(right.coordinates)
      };
    }
  }
  public struct TriVector<V> where V : ISpace {
    public Scalar coordinate;
  }

  #endregion Grassman algebra

  #region Affine space

  public struct Point<V> where V : ISpace {
    public R3Element coordinates;
    // A convex combination of positions is a position.
    public static Point<V> Barycenter(Point<V> q1, Scalar λ1,
                                      Point<V> q2, Scalar λ2) {
      return new Point<V> {
        coordinates = (q1.coordinates * λ1 + q2.coordinates * λ2) / (λ1 + λ2)
      };
    }
  }

  #endregion Affine space

  #region Maps between base spaces

  public struct Rotation<V, W>
    where V : ISpace
    where W : ISpace {
    // The rotation is modeled as a quaternion.
    public Scalar realPart;
    public R3Element imaginaryPart;
    public Rotation<W, V> Inverse {
      get {
        return new Rotation<W, V> {
          realPart = realPart,
          imaginaryPart = -imaginaryPart
        };
      }
    }
    public static R3Element operator *(Rotation<V, W> left,
                                       R3Element right) {
      // TODO(egg): Optimise.
      return Maps.Compose<W, V, W>(
        Maps.Compose<V, V, W>(
          left,
          new Rotation<V, V> { realPart = (Scalar)0, imaginaryPart = right }),
        left.Inverse).imaginaryPart;
    }
    public Vector<W> ActOn(Vector<V> right) {
      return new Vector<W> { coordinates = this * right.coordinates };
    }
  }
  public struct OrthogonalTransformation<V, W>
    where V : ISpace
    where W : ISpace {
    // The orthogonal transformation is modeled as a rotoinversion.
    public Sign determinant;
    public Rotation<V, W> specialOrthogonalMap;

    public static R3Element operator *(OrthogonalTransformation<V, W> left,
                                      R3Element right) {
      return left.specialOrthogonalMap * (left.determinant * right);
    }
    public Vector<W> ActOn(Vector<V> right) {
      return new Vector<W> { coordinates = this * right.coordinates };
    }
    public BiVector<V> ActOn(BiVector<V> right) {
      return new BiVector<V> {
        coordinates = this.specialOrthogonalMap * right.coordinates
      };
    }
    public TriVector<V> ActOn(TriVector<V> right) {
      return new TriVector<V> {
        coordinate = this.determinant * right.coordinate
      };
    }
  }
  public struct EuclideanTransformation<V, W>
    where V : ISpace
    where W : ISpace {
    public OrthogonalTransformation<V, W> orthogonalMap;
    public Vector<W> translation;
  }

  #endregion Maps between base spaces

  #region Methods for the composition of maps

  public static class Maps {
    public static Rotation<V, X> Compose<V, W, X>(Rotation<W, X> left,
                                                  Rotation<V, W> right)
      where V : ISpace
      where W : ISpace
      where X : ISpace {
      return new Rotation<V, X> {
        realPart = left.realPart * right.realPart
                   - left.imaginaryPart.Dot(right.imaginaryPart),
        imaginaryPart = left.realPart * right.imaginaryPart
                        + right.realPart * left.imaginaryPart
                        + left.imaginaryPart.Cross(right.imaginaryPart)
      };
    }
    public static OrthogonalTransformation<V, X> Compose<V, W, X>(
  OrthogonalTransformation<W, X> left,
  OrthogonalTransformation<V, W> right)
      where V : ISpace
      where W : ISpace
      where X : ISpace {
      return new OrthogonalTransformation<V, X> {
        determinant = left.determinant * right.determinant,
        specialOrthogonalMap = Compose<V, W, X>(left.specialOrthogonalMap,
                                                right.specialOrthogonalMap)
      };
    }
    public static EuclideanTransformation<V, X> Compose<V, W, X>(
  EuclideanTransformation<W, X> left,
  EuclideanTransformation<V, W> right)
      where V : ISpace
      where W : ISpace
      where X : ISpace {
      return new EuclideanTransformation<V, X> {
        orthogonalMap = Compose<V, W, X>(left.orthogonalMap,
                                         right.orthogonalMap),
        translation = left.translation
                      + left.orthogonalMap.ActOn(right.translation)
      };
    }
  }

  #endregion Methods for the composition of maps
}