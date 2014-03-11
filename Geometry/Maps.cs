using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

// This namespace is a monstrosity as far as code replication is concerned, but
// there is no way around that short of rewriting the whole thing in C++/CLI.
// This might be done eventually (especially if we need to interface with native
// C++), but is not a high priority. It would make the management of physical
// quantities easier and more systematic.
namespace Geometry {
  // Structs for strong typing in a three dimensional real inner product space
  // and its Grassmann algebra, as well as an affine space.
  // Double-precision floating point numbers, representing reals, are wrapped in
  // the struct Scalar to ensure strong typing. Three-dimensional data is stored
  // as R3Element.
  // We use the Grassman algebra in order to avoid the confusion between pseudo-
  // vectors and vectors and between pseudoscalars and scalars. This allows for
  // handling of indirect changes of coordinates. As a result, there is no cross
  // product except for the underlying data type R3Element.
  // We treat the isomorphism between a space and its dual as implicit, as we
  // are dealing with inner product spaces.
  // TODO(robin): We may want to change this when we get to rotating-pulsating
  // reference frames.
  // A reminder of the equivalents in the Grassman algebra of the cross
  // product; We identify V ^ V with the Lie algebra so(V).
  // . pseudovector = vector x vector, as in L = r x p, becomes the wedge
  //   product,
  //   bivector = vector ^ vector, L = r ^ p.
  // . vector = pseudovector x vector, as in a_Coriolis = 2 Ω x v, becomes
  //   the left action of the bivector,
  //   vector = bivector . vector,  a_Coriolis = 2 Ω . v.
  // . vector = vector x pseudovector, as in F = q v x B, becomes
  //   the dual action,
  //   vector = bivector . vector,  F* = q v* . B.
  // . pseudovector = pseudovector x pseudovector, as in the composition of
  //   rotation vectors Θ12 = α Θ1 + β Θ2 + γ Θ1 x Θ2, becomes the Lie bracket
  //   of so(V).
  //   bivector = [bivector, bivector], Θ12 = α Θ1 + β Θ2 + γ [Θ1, Θ2].
  // . pseudoscalar = (pseudovector|vector), as in Φ = (B|S), becomes the wedge
  //   product,
  //   trivector = bivector ^ vector, Φ = B ^ S.
  // NOTATION:
  // For k-vectors v, w, the inner product (v|w) is [k]vector.InnerProduct(v, w)
  // rather than v.InnerProduct(w) so that the enclosing space is explicit.
  // Compare with the notation v.Wedge(w) for the exterior product of a k-vector
  // v and a l-vector w.
  // We choose to use operators for the vector space operations, thus keeping
  // the enclosing space implicit, i.e., v + w rather than [k]vector.Plus(v, w).
  // The left action a . v of a bivector a on v is a.ActOn(v).
  // The right action v* . a is v.ActedUponBy(a). The commutator is
  // Bivector.Commutator(a, b).

  public interface ISpace { };

  public static class Maps {
    public static Rotation<A, C> Compose<A, B, C>(Rotation<B, C> left,
                                                  Rotation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new Rotation<A, C> {
        RealPart = left.RealPart * right.RealPart
                   - left.ImaginaryPart.Dot(right.ImaginaryPart),
        ImaginaryPart = left.RealPart * right.ImaginaryPart
                        + right.RealPart * left.ImaginaryPart
                        + left.ImaginaryPart.Cross(right.ImaginaryPart)
      };
    }
    public static OrthogonalTransformation<A, C> Compose<A, B, C>(
      OrthogonalTransformation<B, C> left,
      OrthogonalTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new OrthogonalTransformation<A, C> {
        Determinant = left.Determinant * right.Determinant,
        SpecialOrthogonalMap = Compose<A, B, C>(left.SpecialOrthogonalMap,
                                                right.SpecialOrthogonalMap)
      };
    }
    public static EuclideanTransformation<A, C> Compose<A, B, C>(
      EuclideanTransformation<B, C> left,
      EuclideanTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new EuclideanTransformation<A, C> {
        orthogonalMap = Compose<A, B, C>(left.orthogonalMap,
                                         right.orthogonalMap),
        translation = left.translation
                      + left.orthogonalMap.ActOn(right.translation)
      };
    }
    public static RigidTransformation<A, C> Compose<A, B, C>(
      RigidTransformation<B, C> left,
      RigidTransformation<A, B> right)
      where A : ISpace
      where B : ISpace
      where C : ISpace {
      return new RigidTransformation<A, C> {
        orthogonalMap = Compose<A, B, C>(left.orthogonalMap,
                                         right.orthogonalMap),
        translation = left.translation
                      + left.orthogonalMap.ActOn(right.translation)
      };
    }
  }
}