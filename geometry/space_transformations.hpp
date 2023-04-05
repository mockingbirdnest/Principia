#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/conformal_map.hpp"
#include "geometry/orthogonal_map.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace _space_transformations {
namespace internal {

using namespace principia::geometry::_affine_map;
using namespace principia::geometry::_conformal_map;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::quantities::_quantities;

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, Length, OrthogonalMap>;

template<typename FromFrame, typename ToFrame>
using DoubleConformalMap = ConformalMap<double, FromFrame, ToFrame>;

template<typename FromFrame, typename ToFrame>
using Similarity =
    AffineMap<FromFrame, ToFrame, Length, DoubleConformalMap>;

}  // namespace internal

using internal::RigidTransformation;
using internal::Similarity;

}  // namespace _space_transformations
}  // namespace geometry
}  // namespace principia
