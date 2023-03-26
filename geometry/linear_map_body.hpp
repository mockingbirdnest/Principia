#pragma once

#include "geometry/linear_map.hpp"

namespace principia {
namespace geometry {
namespace _linear_map {
namespace internal {

template<typename Map, typename FromFrame, typename ToFrame>
static Map LinearMap<Map, FromFrame, ToFrame>::Identity() {
  return Map::Identity();
}

template<typename Map, typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> LinearMap<Map, FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return static_cast<Map>(*this)(vector);
}

template<typename Map, typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<Map, T>::type
LinearMap<Map, FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<Map, T>::Do(static_cast<Map>(*this), t);
}

template<typename Map, typename FromFrame, typename ToFrame>
void LinearMap<Map, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) {
  FromFrame::WriteToMessage(message->mutable_from_frame());
  ToFrame::WriteToMessage(message->mutable_to_frame());
}

template<typename Map, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
void LinearMap<Map, FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  FromFrame::ReadFromMessage(message.from_frame());
  ToFrame::ReadFromMessage(message.to_frame());
}

}  // namespace internal
}  // namespace _linear_map
}  // namespace geometry
}  // namespace principia
