#pragma once

#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {
namespace _linear_map {
namespace internal {

using namespace principia::base::_mappable;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_grassmann;

template<typename Map, typename FromFrame, typename ToFrame>
class LinearMap {
 public:
  virtual ~LinearMap() = default;

  // The contract of linear maps.  All subclasses must implement these functions
  // lest they go into infinite loops or trigger weird compilation errors.  In
  // C++23, we'll want to use concepts.

  static Map Identity();

  // Our vector are effectively contravariant.  Note that only orthogonal maps
  // operate on bivectors, trivectors, and symmetric bilinear forms because we
  // do not want to distinguish their possible co/contravariance.
  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename T>
  typename Mappable<Map, T>::type operator()(T const& t) const;

  // Apologies for the commented-out code, but we cannot write the return type
  // of this function.
  //
  //   virtual LinearMap<..., ToFrame, FromFrame> Inverse() const = 0;

 protected:
  // Serialization of the frames.  These are just helper functions for
  // implementing the subclasses, they don't dispatch to the subclasses.

  static void WriteToMessage(not_null<serialization::LinearMap*> message);

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<is_serializable_v<F> &&
                                       is_serializable_v<T>>>
  static void ReadFromMessage(serialization::LinearMap const& message);
};

}  // namespace internal

using internal::LinearMap;

}  // namespace _linear_map
}  // namespace geometry
}  // namespace principia

#include "geometry/linear_map_body.hpp"
