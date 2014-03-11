using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
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
}