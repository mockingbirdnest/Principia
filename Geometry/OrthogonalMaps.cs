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
    public static Permutation<A, B> Identity = new Permutation<A, B> {
      BasisImage = CoordinatePermutation.XYZ
    };

    #region ILinearMap implementation

    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B> { Coordinates = this * right.Coordinates };
    }
    public BiVector<B> ActOn(BiVector<A> right) {
      return new BiVector<B> {
        Coordinates = this.Determinant * (this * right.Coordinates)
      };
    }
    public TriVector<B> ActOn(TriVector<A> right) {
      return new TriVector<B> {
        Coordinate = this.Determinant * right.Coordinate
      };
    }

    #endregion ILinearMap implementation

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
  }

  // The rotation is modeled as a quaternion.
  public struct Rotation<A, B> : ILinearMap<A, B>
    where A : ISpace
    where B : ISpace {
    public static Rotation<A, B> Identity = new Rotation<A, B> {
      RealPart = (Scalar)1,
      ImaginaryPart = new R3Element {
        X = (Scalar)0,
        Y = (Scalar)0,
        Z = (Scalar)0
      }
    };

    #region ILinearMap implementation

    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B> { Coordinates = this * right.Coordinates };
    }
    public BiVector<B> ActOn(BiVector<A> right) {
      return new BiVector<B> { Coordinates = this * right.Coordinates };
    }
    public TriVector<B> ActOn(TriVector<A> right) {
      return new TriVector<B> { Coordinate = right.Coordinate };
    }

    #endregion ILinearMap implementation

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
  }

  // The orthogonal transformation is modeled as a rotoinversion.
  public struct OrthogonalTransformation<A, B> : ILinearMap<A, B>
    where A : ISpace
    where B : ISpace {
    public static OrthogonalTransformation<A, B> Identity
      = new OrthogonalTransformation<A, B> {
        Determinant = (Sign)1,
        SpecialOrthogonalMap = Rotation<A, B>.Identity
      };
    public Sign Determinant;
    public Rotation<A, B> SpecialOrthogonalMap;

    public static R3Element operator *(OrthogonalTransformation<A, B> left,
                                      R3Element right) {
      return left.SpecialOrthogonalMap * (left.Determinant * right);
    }

    #region ILinearMap implementation

    public Vector<B> ActOn(Vector<A> right) {
      return new Vector<B> { Coordinates = this * right.Coordinates };
    }
    public BiVector<B> ActOn(BiVector<A> right) {
      return new BiVector<B> {
        Coordinates = this.SpecialOrthogonalMap * right.Coordinates
      };
    }
    public TriVector<B> ActOn(TriVector<A> right) {
      return new TriVector<B> {
        Coordinate = this.Determinant * right.Coordinate
      };
    }

    #endregion ILinearMap implementation
  }
}