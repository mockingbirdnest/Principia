
#pragma once

#include "geometry/identity.hpp"

#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace internal_identity {

template<typename FromFrame, typename ToFrame>
Identity<FromFrame, ToFrame>::Identity() {}

template<typename FromFrame, typename ToFrame>
Sign Identity<FromFrame, ToFrame>::Determinant() const {
  return Sign::Positive();
}

template<typename FromFrame, typename ToFrame>
Identity<ToFrame, FromFrame> Identity<FromFrame, ToFrame>::Inverse() const {
  return Identity<ToFrame, FromFrame>();
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Identity<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>(vector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Identity<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(bivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Identity<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar, template<typename S, typename F> typename Multivector>
SymmetricBilinearForm<Scalar, ToFrame, Multivector>
Identity<FromFrame, ToFrame>::operator()(
    SymmetricBilinearForm<Scalar, FromFrame, Multivector> const& form) const {
  return SymmetricBilinearForm<Scalar, ToFrame, Multivector>(
      form.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<Identity<FromFrame, ToFrame>, T>::type
Identity<FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<Identity, T>::Do(*this, t);
}

template<typename FromFrame, typename ToFrame>
template<template<typename, typename> typename LinearMap>
LinearMap<FromFrame, ToFrame> Identity<FromFrame, ToFrame>::Forget() const {
  return LinearMap<FromFrame, ToFrame>::Identity();
}

template<typename FromFrame, typename ToFrame>
void Identity<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(message->MutableExtension(serialization::Identity::extension));
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Identity<FromFrame, ToFrame> Identity<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Identity::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Identity::extension));
}

template<typename FromFrame, typename ToFrame>
void Identity<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::Identity*> const message) const {}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Identity<FromFrame, ToFrame> Identity<FromFrame, ToFrame>::ReadFromMessage(
    serialization::Identity const& message) {
  return Identity();
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Identity<FromFrame, ToFrame> operator*(
    Identity<ThroughFrame, ToFrame> const& left,
    Identity<FromFrame, ThroughFrame> const& right) {
  return Identity<FromFrame, ToFrame>();
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Identity<FromFrame, ToFrame> const& identity) {
  return out << u8"𝟙";
}

}  // namespace internal_identity
}  // namespace geometry
}  // namespace principia
