using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry {
  public interface ISpace { };

  public static class Maps {
    public static Rotation<A, C> Compose<A, B, C>(Rotation<B, C> left,
                                                  Rotation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new Rotation<A, C>(
         left.RealPart * right.RealPart
                   - left.ImaginaryPart.Dot(right.ImaginaryPart),
         left.RealPart * right.ImaginaryPart
                        + right.RealPart * left.ImaginaryPart
                        + left.ImaginaryPart.Cross(right.ImaginaryPart)
      );
    }
    public static OrthogonalTransformation<A, C> Compose<A, B, C>(
      OrthogonalTransformation<B, C> left,
      OrthogonalTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new OrthogonalTransformation<A, C>(
        left.Determinant * right.Determinant,
        Compose<A, B, C>(left.SpecialOrthogonalMap, right.SpecialOrthogonalMap)
      );
    }
    public static EuclideanTransformation<A, C> Compose<A, B, C>(
      EuclideanTransformation<B, C> left,
      EuclideanTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new EuclideanTransformation<A, C>(
        Compose<A, B, C>(left.OrthogonalMap, right.OrthogonalMap),
        left.Translation + left.OrthogonalMap.ActOn(right.Translation)
      );
    }
    public static RigidTransformation<A, C> Compose<A, B, C>(
      RigidTransformation<B, C> left,
      RigidTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new RigidTransformation<A, C>(
        Compose<A, B, C>(left.SpecialOrthogonalMap, right.SpecialOrthogonalMap),
        left.Translation + left.SpecialOrthogonalMap.ActOn(right.Translation)
      );
    }
  }
}