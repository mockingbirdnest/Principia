#pragma once

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> LinearMap>
class AffineMap {
 public:
  typedef Vector<Scalar, FromFrame> FromVector;
  typedef Vector<Scalar, ToFrame> ToVector;
  AffineMap(Point<FromVector> from_origin,
            Point<ToVector> to_origin,
            LinearMap<FromFrame, ToFrame> linear_map);

  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> Inverse() const;
  Point<ToVector> operator()(Point<FromVector> point) const;
 private:
  ToVector translation_;
  LinearMap<FromFrame, ToFrame> linear_map_;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& right);

}  // namespace geometry
}  // namespace principia
