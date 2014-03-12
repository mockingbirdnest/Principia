using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  public struct Point<A> where A : ISpace {
    public readonly R3Element Coordinates;
    public Point(R3Element coordinates) {
      Coordinates = coordinates;
    }
    // A convex combination of positions is a position.
    public static Point<A> Barycenter(Point<A> q1, Scalar λ1,
                                      Point<A> q2, Scalar λ2) {
      return new Point<A>((q1.Coordinates * λ1 + q2.Coordinates * λ2) / (λ1 + λ2));
    }
  }
}