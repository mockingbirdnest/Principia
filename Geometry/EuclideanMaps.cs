using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
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
}