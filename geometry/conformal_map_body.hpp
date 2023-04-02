#pragma once

#include "geometry/conformal_map.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _conformal_map {
namespace internal {

using namespace principia::quantities::_elementary_functions;

template<typename Scalar, typename FromFrame, typename ToFrame>
Cube<Scalar> ConformalMap<Scalar, FromFrame, ToFrame>::Determinant() const {
  return Pow<3>(scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
ConformalMap<Inverse<Scalar>, ToFrame, FromFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::Inverse() const {
  return ConformalMap<quantities::Inverse<Scalar>, ToFrame, FromFrame>(
    1 / scale_, quaternion_.Inverse());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::Identity() {
  return ConformalMap(1, Quaternion(1));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename VScalar>
Vector<Product<VScalar, Scalar>, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::operator()(
    Vector<VScalar, FromFrame> const& vector) const {
  return MakeHomothecy()(MakeRotation()(MakeSignature()(vector)));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename T>
typename Mappable<ConformalMap<Scalar, FromFrame, ToFrame>, T>::type
ConformalMap<Scalar, FromFrame, ToFrame>::operator()(T const& t) const {
  return Mappable<ConformalMap, T>::Do(*this, t);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void ConformalMap<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<ConformalMap, FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::ConformalMap::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<ConformalMap, FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::ConformalMap::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::ConformalMap::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void ConformalMap<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::ConformalMap*> const message) const {
  scale_.WriteToMessage(message->mutable_scale());
  quaternion_.WriteToMessage(message->mutable_quaternion());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
ConformalMap<Scalar, FromFrame, ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::ConformalMap const& message) {
  return ConformalMap(Scalar::ReadFromMessage(message.scale()),
                      Quaternion::ReadFromMessage(message.quaternion()));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
ConformalMap<Scalar, FromFrame, ToFrame>::ConformalMap(
    Scalar const& scale,
    Quaternion const& quaternion)
    : scale_(scale),
      quaternion_(quaternion) {
  CHECK_LT(Scalar{}, scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
constexpr auto ConformalMap<Scalar, FromFrame, ToFrame>::MakeSignature() ->
Signature<FromFrame, SignedFrame> {
  if constexpr (FromFrame::handedness == ToFrame::handedness) {
    return Signature<FromFrame, SignedFrame>::Identity();
  } else {
    return Signature<FromFrame, SignedFrame>::CentralInversion();
  }
}

template<typename Scalar, typename FromFrame, typename ToFrame>
auto ConformalMap<Scalar, FromFrame, ToFrame>::MakeRotation() const ->
Rotation<SignedFrame, RotatedAndSignedFrame> {
  return Rotation<SignedFrame, RotatedAndSignedFrame>(quaternion_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
auto ConformalMap<Scalar, FromFrame, ToFrame>::MakeHomothecy() const ->
Homothecy<Scalar, RotatedAndSignedFrame, ToFrame> {
  using H = Homothecy<Scalar, RotatedAndSignedFrame, ToFrame>;
  return H(typename H::PrivateConstructor{}, scale_);
}

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
ConformalMap<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    ConformalMap<LScalar, ThroughFrame, ToFrame> const& left,
    ConformalMap<RScalar, FromFrame, ThroughFrame> const& right) {
  using Map = ConformalMap<Product<LScalar, RScalar>, FromFrame, ToFrame>;
  return Map(
      left.scale_ * right.scale_,
      left.quaternion_ * right.quaternion_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    ConformalMap<Scalar, FromFrame, ToFrame> const& conformal_map) {
  return out << "{scale: " << conformal_map.scale_
             << ", quaternion: " << conformal_map.quaternion_ << "}";
}

}  // namespace internal
}  // namespace _conformal_map
}  // namespace geometry
}  // namespace principia
