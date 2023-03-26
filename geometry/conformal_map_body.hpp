#pragma once

#include "geometry/conformal_map.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _conformal_map {
namespace internal {

using namespace principia::quantities::_elementary_functions;

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>::ConformalMap(
    Scalar const& scale,
    OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map)
    : ConformalMap(PrivateConstructor{}, scale, orthogonal_map) {}

template<typename Scalar, typename FromFrame, typename ToFrame>
Cube<Scalar> ConformalMap<Scalar, FromFrame, ToFrame>::Determinant() const {
  return Pow<3>(scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
ConformalMap<Inverse<Scalar>, ToFrame, FromFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::Inverse() const {
  return ConformalMap<quantities::Inverse<Scalar>, ToFrame, FromFrame>(
    1 / scale_, orthogonal_map_.Inverse());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::Identity() {
  return ConformalMap(
      PrivateConstructor{}, 1, OrthogonalMap<FromFrame, ToFrame>::Identity());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename VScalar>
Vector<Product<VScalar, Scalar>, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::operator()(
    Vector<VScalar, FromFrame> const& vector) const {
  return Vector<Product<VScalar, Scalar>, ToFrame>(
      orthogonal_map_(vector).coordinates() * scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<ConformalMap<Scalar, FromFrame, ToFrame>, T>::type
ConformalMap<Scalar, FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<ConformalMap, T>::Do(*this, t);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void ConformalMap<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::ConformalMap::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::ConformalMap::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::ConformalMap::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void ConformalMap<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::ConformalMap*> const message) const {
  scale_.WriteToMessage(message->mutable_scale());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::ConformalMap const& message) {
  return ConformalMap(Scalar::ReadFromMessage(message.scale()));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::ConformalMap(
    PrivateConstructor,
    Scalar const& scale,
    OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map)
    : scale_(scale),
      orthogonal_map_(orthogonal_map) {
  CHECK_LT(Scalar{}, scale_);
}

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
ConformalMap<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    ConformalMap<LScalar, ThroughFrame, ToFrame> const& left,
    ConformalMap<RScalar, FromFrame, ThroughFrame> const& right) {
  return ConformalMap<Product<LScalar, RScalar>, FromFrame, ToFrame>(
      left.scale_ * right.scale_,
      left.orthogonal_map_ * right.orthogonal_map_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    ConformalMap<Scalar, FromFrame, ToFrame> const& conformal_map) {
  return out << "{scale: " << conformal_map.scale_
             << ", orthogonal_map: " << conformal_map.orthogonal_map_ << "}";
}

}  // namespace internal
}  // namespace _conformal_map
}  // namespace geometry
}  // namespace principia
