
#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_sphere {

using quantities::Length;
using quantities::Square;

template<typename Frame>
class Sphere {
 public:
  Sphere(Position<Frame> const& centre,
         Length const& radius);

  Position<Frame> const& centre() const;
  Length const& radius() const;
  Square<Length> const& radius²() const;

 private:
  Position<Frame> const centre_;
  Length const radius_;
  Square<Length> const radius²_;
};

}  // namespace internal_sphere

using internal_sphere::Sphere;

}  // namespace geometry
}  // namespace principia

#include "geometry/sphere_body.hpp"
