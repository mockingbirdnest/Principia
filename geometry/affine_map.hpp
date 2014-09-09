#pragma once

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
class AffineMap {
 public:
  using FromVector = Vector<Scalar, FromFrame>;
  using ToVector = Vector<Scalar, ToFrame>;
  // The identity map.
  AffineMap() = default;
  AffineMap(Point<FromVector> const& from_origin,
            Point<ToVector> const& to_origin,
            LinearMap<FromFrame, ToFrame> const& linear_map);

  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> Inverse() const;
  Point<ToVector> operator()(Point<FromVector> const& point) const;

 private:
  // The map is internally represented as x ↦ linear_map_(x) + translation_.
  LinearMap<FromFrame, ToFrame> linear_map_;
  ToVector translation_;

  template<typename From, typename Through, typename To, typename S,
           template<typename, typename> class Map>
  friend AffineMap<From, To, S, Map> operator*(
      AffineMap<Through, To, S, Map> const& left,
      AffineMap<From, To, S, Map> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/affine_map_body.hpp"
