#pragma once

#include "geometry/homothecy.hpp"

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace _homothecy {
namespace internal {

using namespace principia::quantities::_elementary_functions;

template<typename Scalar, typename FromFrame, typename ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::Homothecy(Scalar const& scale)
    : scale_(scale) {}

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
template<typename VScalar>
Vector<Product<VScalar, Scalar>, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::operator()(
    Vector<VScalar, FromFrame> const& vector) const {
  return Vector<Product<VScalar, Scalar>, ToFrame>(vector.coordinates() *
                                                   scale_);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename BScalar>
Bivector<Product<BScalar, Square<Scalar>>, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::operator()(
    Bivector<BScalar, FromFrame> const& bivector) const {
  return Bivector<Product<BScalar, Scalar>, ToFrame>(bivector.coordinates() *
                                                     Pow<2>(scale_));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename TScalar>
Trivector<Product<TScalar, Cube<Scalar>>, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::operator()(
    Trivector<TScalar, FromFrame> const& trivector) const {
  return Trivector<Product<TScalar, Cube<Scalar>>, ToFrame>(
      trivector.coordinates() * Pow<3>(scale_));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename Scalar, typename Scalar, template<typename Scalar, typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, ToFrame, Multivector>
Homothecy<Scalar, FromFrame, ToFrame>::operator()(
    SymmetricBilinearForm<Scalar, FromFrame, Multivector> const& form) const {
  return SymmetricBilinearForm<Scalar, ToFrame, Multivector>(
      form.coordinates());
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<Homothecy<Scalar, FromFrame, ToFrame>, T>::type
Homothecy<Scalar, FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<Homothecy, T>::Do(*this, t);
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<template<typename, typename> typename LinearMap>
LinearMap<FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::Forget() const {
  return LinearMap<FromFrame, ToFrame>::Homothecy();
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void Homothecy<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::Homothecy::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Homothecy<Scalar, FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Homothecy::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Homothecy::extension));
}

template<typename Scalar, typename FromFrame, typename ToFrame>
void Homothecy<Scalar, FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::Homothecy*> const message) const {}

template<typename Scalar, typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Homothecy<Scalar, FromFrame, ToFrame>
Homothecy<Scalar, FromFrame, ToFrame>::ReadFromMessage(
    serialization::Homothecy const& message) {
  return Homothecy();
}

template<typename Scalar, typename FromFrame, typename ThroughFrame, typename ToFrame>
Homothecy<Scalar, FromFrame, ToFrame> operator*(
    Homothecy<Scalar, ThroughFrame, ToFrame> const& left,
    Homothecy<Scalar, FromFrame, ThroughFrame> const& right) {
  return Homothecy<Scalar, FromFrame, ToFrame>();
}

template<typename Scalar, typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Homothecy<Scalar, FromFrame, ToFrame> const& homothecy) {
  return out << "𝟙";
}

}  // namespace internal
}  // namespace _homothecy
}  // namespace geometry
}  // namespace principia
