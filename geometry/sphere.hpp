
#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_sphere {

using quantities::Product;

template<typename Scalar, typename Frame>
class Sphere {
 public:
  Sphere(Point<Vector<Scalar, Frame>> const& centre,
         Scalar const& radius);

  Point<Vector<Scalar, Frame>> const& centre() const;
  Scalar const& radius() const;
  Product<Scalar, Scalar> const radius²() const;

 private:
  Point<Vector<Scalar, Frame>> const centre_;
  Scalar const radius_;
  Product<Scalar, Scalar> const radius²_;
};

}  // namespace internal_sphere

using internal_sphere::Sphere;

}  // namespace geometry
}  // namespace principia

#include "geometry/sphere_body.hpp"
