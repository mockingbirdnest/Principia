using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  internal interface IAffineMap<A, B>
    where A : ISpace
    where B : ISpace {
    Point<B> ActOn(Point<A> right);
  }
  public struct RigidTransformation<A, B> : IAffineMap<A, B>
    where A : ISpace
    where B : ISpace {
    public Rotation<A, B> SpecialOrthogonalMap;
    public Vector<B> Translation;

    #region IAffineMap implementation

    public Point<B> ActOn(Point<A> right) {
      return Translation.Translate(
        new Point<B> {
          Coordinates = SpecialOrthogonalMap * right.Coordinates
        });
    }

    #endregion IAffineMap implementation

    public EuclideanTransformation<A, B> Forget() {
      return new EuclideanTransformation<A, B> {
        OrthogonalMap = SpecialOrthogonalMap.Forget(),
        Translation = Translation
      };
    }
  }
  public struct EuclideanTransformation<A, B> : IAffineMap<A, B>
    where A : ISpace
    where B : ISpace {
    public OrthogonalTransformation<A, B> OrthogonalMap;
    public Vector<B> Translation;

    #region IAffineMap implementation

    public Point<B> ActOn(Point<A> right) {
      return Translation.Translate(
        new Point<B> {
          Coordinates = OrthogonalMap * right.Coordinates
        });
    }

    #endregion IAffineMap implementation
  }
}