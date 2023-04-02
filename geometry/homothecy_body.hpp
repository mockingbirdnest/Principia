#pragma once

#include "geometry/homothecy.hpp"

#include "geometry/quaternion.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _homothecy {
namespace internal {

using namespace principia::geometry::_quaternion;
using namespace principia::quantities::_elementary_functions;

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename>
Homothecy<Scalar, FromFrame, ToFrame>::Homothecy(Scalar const& scale)
    : Homothecy(PrivateConstructor{}, scale) {}

template<typename Scalar, typename FromFrame, typename ToFrame>
Cube<Scalar> Homothecy<Scalar, FromFrame, ToFrame>::Determinant() const {
  return Pow<3>(scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
Homothecy<Inverse<Scalar>, ToFrame, FromFrame>
Homothecy<Scalar, FromFrame, ToFrame>::Inverse() const {
  return Homothecy<quantities::Inverse<Scalar>, ToFrame, FromFrame>(
    1 / scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename>
Homothecy<Scalar, FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::Identity() {
  return Homothecy(PrivateConstructor{}, 1);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename VScalar>
Vector<Product<VScalar, Scalar>, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::operator()(
    Vector<VScalar, FromFrame> const& vector) const {
  return Vector<Product<VScalar, Scalar>, ToFrame>(vector.coordinates() *
                                                   scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename T>
typename Mappable<Homothecy<Scalar, FromFrame, ToFrame>, T>::type
Homothecy<Scalar, FromFrame, ToFrame>::operator()(T const& t) const {
  return Mappable<Homothecy, T>::Do(*this, t);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<template<typename, typename, typename> typename ConformalMap>
ConformalMap<Scalar, FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::Forget() const {
  return ConformalMap<Scalar, FromFrame, ToFrame>(scale_, Quaternion(1));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void Homothecy<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<Homothecy, FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::Homothecy::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Homothecy<Scalar, FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<Homothecy, FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Homothecy::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Homothecy::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void Homothecy<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::Homothecy*> const message) const {
  scale_.WriteToMessage(message->mutable_scale());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Homothecy<Scalar, FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::Homothecy const& message) {
  return Homothecy(PrivateConstructor{},
                   Scalar::ReadFromMessage(message.scale()));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::Homothecy(PrivateConstructor,
                                                 Scalar const& scale)
    : scale_(scale) {
  CHECK_LT(Scalar{}, scale_);
}

template<typename LScalar, typename RScalar,
         typename FromFrame, typename ThroughFrame, typename ToFrame>
Homothecy<Product<LScalar, RScalar>, FromFrame, ToFrame> operator*(
    Homothecy<LScalar, ThroughFrame, ToFrame> const& left,
    Homothecy<RScalar, FromFrame, ThroughFrame> const& right) {
  return Homothecy<Product<LScalar, RScalar>, FromFrame, ToFrame>(left.scale_ *
                                                                  right.scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    Homothecy<Scalar, FromFrame, ToFrame> const& homothecy) {
  return out << homothecy.scale_;
}

}  // namespace internal
}  // namespace _homothecy
}  // namespace geometry
}  // namespace principia
