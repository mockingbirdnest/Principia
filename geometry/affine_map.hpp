#pragma once

#include "base/concepts.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _affine_map {
namespace internal {

using namespace principia::base::_concepts;
using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_point;

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap_>
class AffineMap final {
 public:
  template<typename F, typename T>
  using LinearMap = LinearMap_<F, T>;
  using FromVector = Vector<Scalar, FromFrame>;
  using ToVector = Vector<Scalar, ToFrame>;

  AffineMap(Point<FromVector> const& from_origin,
            Point<ToVector> const& to_origin,
            LinearMap<FromFrame, ToFrame> linear_map);

  AffineMap<ToFrame, FromFrame, Scalar, LinearMap_> Inverse() const;
  Point<ToVector> operator()(Point<FromVector> const& point) const;

  template<template<typename, typename> typename OtherAffineMap>
  OtherAffineMap<FromFrame, ToFrame> Forget() const;

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<F::handedness == T::handedness>>
  static AffineMap Identity();

  LinearMap<FromFrame, ToFrame> const& linear_map() const;

  void WriteToMessage(not_null<serialization::AffineMap*> message) const;
  static AffineMap ReadFromMessage(serialization::AffineMap const& message)
    requires serializable<FromFrame> && serializable<ToFrame>;

 private:
  Point<FromVector> from_origin_;
  Point<ToVector> to_origin_;
  LinearMap<FromFrame, ToFrame> linear_map_;

  template<typename From, typename Through, typename To, typename S,
           template<typename, typename> class Map>
  friend AffineMap<From, To, S, Map> operator*(
      AffineMap<Through, To, S, Map> const& left,
      AffineMap<From, Through, S, Map> const& right);
  template<typename From, typename To, typename S,
           template<typename, typename> class Map>
  friend std::ostream& operator<<(
      std::ostream& out,
      AffineMap<From, To, S, Map> const& affine_map);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ThroughFrame, Scalar, LinearMap> const& right);

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
std::ostream& operator<<(
    std::ostream& out,
    AffineMap<FromFrame, ToFrame, Scalar, LinearMap> const& affine_map);

}  // namespace internal

using internal::AffineMap;

}  // namespace _affine_map
}  // namespace geometry
}  // namespace principia

#include "geometry/affine_map_body.hpp"
