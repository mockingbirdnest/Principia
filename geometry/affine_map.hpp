#pragma once

#include "geometry/point.hpp"
#include "geometry/grassmann.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
class AffineMap {
 public:
  using FromVector = Vector<Scalar, FromFrame>;
  using ToVector = Vector<Scalar, ToFrame>;

  AffineMap(Point<FromVector> const& from_origin,
            Point<ToVector> const& to_origin,
            LinearMap<FromFrame, ToFrame> const& linear_map);

  AffineMap(AffineMap const&) = default;

  AffineMap<ToFrame, FromFrame, Scalar, LinearMap> Inverse() const;
  Point<ToVector> operator()(Point<FromVector> const& point) const;

  static AffineMap Identity();

  LinearMap<FromFrame, ToFrame> const& linear_map() const;

  void WriteToMessage(not_null<serialization::AffineMap*> const message) const;
  static AffineMap ReadFromMessage(serialization::AffineMap const& message);

 private:
  AffineMap() = default;

  Point<FromVector> from_origin_;
  Point<ToVector> to_origin_;
  LinearMap<FromFrame, ToFrame> linear_map_;

  template<typename From, typename Through, typename To, typename S,
           template<typename, typename> class Map>
  friend AffineMap<From, To, S, Map> operator*(
      AffineMap<Through, To, S, Map> const& left,
      AffineMap<From, Through, S, Map> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ThroughFrame, Scalar, LinearMap> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/affine_map_body.hpp"
