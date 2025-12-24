#pragma once

#include "base/algebra.hpp"
#include "geometry/space.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _sphere {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::geometry::_space;
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

#include "geometry/sphere_body.hpp"
