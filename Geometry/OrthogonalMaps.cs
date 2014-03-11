using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
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
}