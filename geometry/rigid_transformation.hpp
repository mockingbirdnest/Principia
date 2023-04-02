#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/orthogonal_map.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _rigid_transformation {
namespace internal {

using namespace principia::geometry::_affine_map;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::quantities::_quantities;

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, Length, OrthogonalMap>;

}  // namespace internal

using internal::RigidTransformation;

}  // namespace _rigid_transformation
}  // namespace geometry
}  // namespace principia
