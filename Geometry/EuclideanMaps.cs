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
    public readonly Rotation<A, B> SpecialOrthogonalMap;
    public readonly Vector<B> Translation;

    public RigidTransformation(Rotation<A, B> specialOrthogonalMap,
                               Vector<B> translation) {
      SpecialOrthogonalMap = specialOrthogonalMap;
      Translation = translation;
    }

    #region IAffineMap implementation

    public Point<B> ActOn(Point<A> right) {
      return Translation.Translate(new Point<B>(SpecialOrthogonalMap *
                                                right.Coordinates));
    }

    #endregion IAffineMap implementation

    public EuclideanTransformation<A, B> Forget() {
      return new EuclideanTransformation<A, B>(SpecialOrthogonalMap.Forget(),
                                               Translation);
    }
  }
  public struct EuclideanTransformation<A, B> : IAffineMap<A, B>
    where A : ISpace
    where B : ISpace {
    public readonly OrthogonalTransformation<A, B> OrthogonalMap;
    public readonly Vector<B> Translation;
    public EuclideanTransformation(OrthogonalTransformation<A, B> orthogonalMap,
                                   Vector<B> translation) {
      OrthogonalMap = orthogonalMap;
      Translation = translation;
    }

    #region IAffineMap implementation

    public Point<B> ActOn(Point<A> right) {
      return Translation.Translate(new Point<B>(OrthogonalMap *
                                                right.Coordinates));
    }

    #endregion IAffineMap implementation
  }
}