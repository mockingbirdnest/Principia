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
    private readonly Vector<B> translation;
    private readonly Rotation<A, B> specialOrthogonalMap;
    public RigidTransformation(Rotation<A, B> specialOrthogonalMap,
                               Vector<B> translation) {
      this.specialOrthogonalMap = specialOrthogonalMap;
      this.translation = translation;
    }
    public Rotation<A, B> SpecialOrthogonalMap {
      get { return specialOrthogonalMap; }
    }
    public Vector<B> Translation { get { return translation; } }

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
    private readonly OrthogonalTransformation<A, B> orthogonalMap;
    private readonly Vector<B> translation;
    public EuclideanTransformation(OrthogonalTransformation<A, B> orthogonalMap,
                                   Vector<B> translation) {
      this.orthogonalMap = orthogonalMap;
      this.translation = translation;
    }
    public OrthogonalTransformation<A, B> OrthogonalMap {
      get { return orthogonalMap; }
    }
    public Vector<B> Translation { get { return translation; } }

    #region IAffineMap implementation

    public Point<B> ActOn(Point<A> right) {
      return Translation.Translate(new Point<B>(OrthogonalMap *
                                                right.Coordinates));
    }

    #endregion IAffineMap implementation
  }
}