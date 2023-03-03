#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace geometry {
namespace _sphere {
namespace internal {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;

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

}  // namespace internal

using internal::Sphere;

}  // namespace _sphere
}  // namespace geometry
}  // namespace principia

namespace principia::geometry {
using namespace principia::geometry::_sphere;
}  // namespace principia::geometry

#include "geometry/sphere_body.hpp"
