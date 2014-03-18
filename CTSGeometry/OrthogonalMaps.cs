using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  // The linear maps between base spaces induce unique maps between Grassmann
  // algebras.
  public interface ILinearMap<A, B>
    where A : ISpace
    where B : ISpace {
    Vector<B> ActOn(Vector<A> right);
    BiVector<B> ActOn(BiVector<A> right);
    TriVector<B> ActOn(TriVector<A> right);
  }
  // A permutation of the coordinates. Obviously not coordinate-free, but
  // practical. There are no precision losses when composing or applying
  // permutations.
  public struct Permutation<A, B> : ILinearMap<A, B>
    where A : ISpace
    where B : ISpace {
    // Odd has the sign bit set to 1 and all other bits set to 0.
    private const int x = 0, y = 1, z = 2, index = 3, odd = Int32.MinValue;
    private static readonly Scalar SqrtHalf = Scalar.Sqrt((Scalar).5);
    private static readonly R3Element[] quaternionImaginaryParts
      = new R3Element[] {
        new R3Element((Scalar)0,(Scalar)0,(Scalar)0),
        new R3Element((Scalar).5, (Scalar).5, (Scalar).5),
        new R3Element((Scalar).5, (Scalar).5, (Scalar).5),
        new R3Element((Scalar)0, -SqrtHalf, SqrtHalf),
        new R3Element(-SqrtHalf, SqrtHalf, (Scalar)0),
        new R3Element(-SqrtHalf, (Scalar)0, SqrtHalf)
      };
    private static readonly Scalar[] quaternionRealParts = new Scalar[] {
      (Scalar)1, (Scalar)(-.5), (Scalar).5, (Scalar)0, (Scalar)0, (Scalar)0
    };
    private readonly CoordinatePermutation basisImage;
    public Permutation(CoordinatePermutation basisImage) {
      this.basisImage = basisImage;
    }
    public enum CoordinatePermutation {
      XYZ = (x << x * 2) + (y << y * 2) + (z << z * 2) + (0 << index * 2),
      YZX = (y << x * 2) + (z << y * 2) + (x << z * 2) + (1 << index * 2),
      ZXY = (z << x * 2) + (x << y * 2) + (y << z * 2) + (2 << index * 2),
      XZY = (x << x * 2) + (z << y * 2) + (y << z * 2) + odd + (3 << index * 2),
      ZYX = (z << x * 2) + (y << y * 2) + (x << z * 2) + odd + (4 << index * 2),
      YXZ = (y << x * 2) + (x << y * 2) + (z << z * 2) + odd + (5 << index * 2)
    }
    public static Permutation<A, B> Identity {
      get { return new Permutation<A, B>(CoordinatePermutation.XYZ); }
    }

    #region ILinearMap implementation

    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B>(this * right.Coordinates);
    }
    public BiVector<B> ActOn(BiVector<A> right) {
      return new BiVector<B>(this.Determinant * (this * right.Coordinates));
    }
    public TriVector<B> ActOn(TriVector<A> right) {
      return new TriVector<B>(this.Determinant * right.Coordinate);
    }

    #endregion ILinearMap implementation

    public Sign Determinant {
      get { return (Sign)(int)basisImage; }
    }
    public static R3Element operator *(Permutation<A, B> left,
                                       R3Element right) {
      R3Element result = new R3Element();
      for (int i = 0; i < 3; ++i) {
        result[i] = right[((3 << i * 2) | (int)left.basisImage) >> i * 2];
      }
      return result;
    }
    public OrthogonalTransformation<A, B> Forget() {
      int i = ((8 << index * 2) | (int)basisImage) >> index * 2;
      return new OrthogonalTransformation<A, B>(
        this.Determinant,
        new Rotation<A, B>(quaternionRealParts[i], quaternionImaginaryParts[i])
        );
    }
  }

  // The rotation is modeled as a quaternion.
  public struct Rotation<A, B> : ILinearMap<A, B>
    where A : ISpace
    where B : ISpace {
    private readonly R3Element imaginaryPart;
    private readonly Scalar realPart;
    public Rotation(Scalar realPart, R3Element imaginaryPart) {
      this.realPart = realPart;
      this.imaginaryPart = imaginaryPart;
    }
    public static Rotation<A, B> Identity {
      get {
        return new Rotation<A, B>((Scalar)1,
                                  new R3Element((Scalar)0,
                                                (Scalar)0,
                                                (Scalar)0));
      }
    }
    public Scalar RealPart { get { return realPart; } }
    public R3Element ImaginaryPart { get { return imaginaryPart; } }
    public static R3Element operator *(Rotation<A, B> left,
                                       R3Element right) {
      // TODO(egg): Optimise.
      return Maps.Compose<B, A, B>(
        Maps.Compose<A, A, B>(left, new Rotation<A, A>((Scalar)0, right)),
        left.Inverse()).ImaginaryPart;
    }
    public Rotation<B, A> Inverse() {
      return new Rotation<B, A>(RealPart, -ImaginaryPart);
    }

    #region ILinearMap implementation

    public BiVector<B> ActOn(BiVector<A> right) {
      return new BiVector<B>(this * right.Coordinates);
    }

    public TriVector<B> ActOn(TriVector<A> right) {
      return new TriVector<B>(right.Coordinate);
    }
    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B>(this * right.Coordinates);
    }

    #endregion ILinearMap implementation

    public OrthogonalTransformation<A, B> Forget() {
      return new OrthogonalTransformation<A, B>((Sign)1, this);
    }
  }

  // The orthogonal transformation is modeled as a rotoinversion.
  public struct OrthogonalTransformation<A, B> : ILinearMap<A, B>
    where A : ISpace
    where B : ISpace {
    private readonly Sign determinant;
    private readonly Rotation<A, B> specialOrthogonalMap;
    public OrthogonalTransformation(Sign determinant, Rotation<A, B> specialOrthogonalMap) {
      this.determinant = determinant;
      this.specialOrthogonalMap = specialOrthogonalMap;
    }
    public static OrthogonalTransformation<A, B> Identity {
      get {
        return new OrthogonalTransformation<A, B>((Sign)1,
                                                  Rotation<A, B>.Identity);
      }
    }
    public Sign Determinant { get { return determinant; } }
    public Rotation<A, B> SpecialOrthogonalMap {
      get { return specialOrthogonalMap; }
    }
    public static R3Element operator *(OrthogonalTransformation<A, B> left,
                                      R3Element right) {
      return left.SpecialOrthogonalMap * (left.Determinant * right);
    }

    #region ILinearMap implementation

    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B>(this * right.Coordinates);
    }
    public BiVector<B> ActOn(BiVector<A> right) {
      return new BiVector<B>(this.SpecialOrthogonalMap * right.Coordinates);
    }
    public TriVector<B> ActOn(TriVector<A> right) {
      return new TriVector<B>(this.Determinant * right.Coordinate);
    }

    #endregion ILinearMap implementation
  }
}