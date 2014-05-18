#pragma once

#include "geometry/affine_space.hpp"

namespace principia {
namespace geometry {

template<typename FromVector, typename ToVector, typename LinearMap>
class AffineMap {
 public:
  AffineMap(Point<FromVector> from_origin,
            Point<ToVector> to_origin,
            LinearMap linear_map);
  Point<ToVector> operator()(Point<FromVector> point);
 private:
  ToVector translation_;
  LinearMap linear_map_;
};

}  // namespace geometry
}  // namespace principia
