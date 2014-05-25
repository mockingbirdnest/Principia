#pragma once

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
class AffineMap {
 public:
  typedef Vector<Scalar, FromFrame> FromVector;
  typedef Vector<Scalar, ToFrame> ToVector;
  // The identity map.
  AffineMap() = default;
  AffineMap(Point<FromVector> const& from_origin,
            Point<ToVector> const& to_origin,
            LinearMap<FromFrame, ToFrame> const& linear_map);

  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> Inverse() const;
  Point<ToVector> operator()(Point<FromVector> const& point) const;

 private:
  // The map is internally represented as x -> linear_map_(x) + translation_.
  ToVector translation_;
  LinearMap<FromFrame, ToFrame> linear_map_;
  
  template<typename FromFrame, typename ThroughFrame, typename ToFrame,
           typename Scalar, template<typename, typename> class LinearMap>
  friend AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
      AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
      AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/affine_map_body.hpp"
