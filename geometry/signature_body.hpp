#pragma once

#include "geometry/signature.hpp"

#include "base/traits.hpp"

namespace principia {
namespace geometry {
namespace internal_signature {

using base::is_same_template_v;

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   Sign const y,
                                                   DeduceSign const z)
    : x_(x), y_(y), z_(x * y * determinant_) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   DeduceSign const y,
                                                   Sign const z)
    : x_(x), y_(x * determinant_ * z), z_(z) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(DeduceSign const x,
                                                   Sign const y,
                                                   Sign const z)
    : x_(determinant_ * y * z), y_(y), z_(z) {}

template<typename FromFrame, typename ToFrame>
constexpr Sign Signature<FromFrame, ToFrame>::Determinant() const {
  return determinant_;
}

template<typename FromFrame, typename ToFrame>
constexpr Signature<ToFrame, FromFrame> Signature<FromFrame, ToFrame>::Inverse()
    const {
  return Signature<ToFrame, FromFrame>(x_, y_, z_);
}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
constexpr Signature<FromFrame, ToFrame>
Signature<FromFrame, ToFrame>::Identity() {
  return Signature(
      Sign::Positive(), Sign::Positive(), DeduceSignPreservingOrientation{});
}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
constexpr Signature<FromFrame, ToFrame>
Signature<FromFrame, ToFrame>::CentralInversion() {
  return Signature(
      Sign::Negative(), Sign::Negative(), DeduceSignReversingOrientation{});
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Signature<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>({x_ * vector.coordinates().x,
                                  y_ * vector.coordinates().y,
                                  z_ * vector.coordinates().z});
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Signature<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(
      {determinant_ * x_ * bivector.coordinates().x,
       determinant_ * y_ * bivector.coordinates().y,
       determinant_ * z_ * bivector.coordinates().z});
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Signature<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(determinant_ * trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar, template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, ToFrame, Multivector>
Signature<FromFrame, ToFrame>::operator()(
    SymmetricBilinearForm<Scalar, FromFrame, Multivector> const& form) const {
  // TODO(egg): This should be writeable as "*= x_ * y_".
  auto coordinates = form.coordinates();
  if (x_ != y_) {
    coordinates(0, 1) *= -1;
    coordinates(1, 0) *= -1;
  }
  if (y_ != z_) {
    coordinates(1, 2) *= -1;
    coordinates(2, 1) *= -1;
  }
  if (z_ != x_) {
    coordinates(2, 0) *= -1;
    coordinates(0, 2) *= -1;
  }
  return SymmetricBilinearForm<Scalar, ToFrame, Multivector>(coordinates);
}

template<typename FromFrame, typename ToFrame>
template<template<typename, typename> typename LinearMap>
LinearMap<FromFrame, ToFrame> Signature<FromFrame, ToFrame>::Forget() const {
  static_assert(is_same_template_v<LinearMap, OrthogonalMap> ||
                is_same_template_v<LinearMap, Rotation>,
                "Unable to forget signature");
  Quaternion quaternion;
  if (x_ == y_ && y_ == z_) {
    quaternion = Quaternion(1);
  } else {
    // The signature is neither +++ nor ---, so dividing it by its determinant
    // yields a 180° rotation around one of the axes (+--, -+-, or --+).
    R3Element<double> const axis = (determinant_ * x_).is_positive()
                                        ? R3Element<double>{1, 0, 0}
                                        : (determinant_ * y_).is_positive()
                                              ? R3Element<double>{0, 1, 0}
                                              : R3Element<double>{0, 0, 1};
    quaternion = Quaternion(0, axis);
  }
  return LinearMap<FromFrame, ToFrame>(quaternion);
}

template<typename FromFrame, typename ToFrame>
void Signature<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::Signature::extension));
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Signature<FromFrame, ToFrame> Signature<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Signature::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Signature::extension));
}

template<typename FromFrame, typename ToFrame>
void Signature<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::Signature*> const message) const {
  x_.WriteToMessage(message->mutable_x());
  y_.WriteToMessage(message->mutable_y());
  z_.WriteToMessage(message->mutable_z());
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Signature<FromFrame, ToFrame> Signature<FromFrame, ToFrame>::ReadFromMessage(
    serialization::Signature const& message) {
  auto const x = Sign::ReadFromMessage(message.x());
  auto const y = Sign::ReadFromMessage(message.y());
  auto const z = Sign::ReadFromMessage(message.z());
  CHECK_EQ(x * y * z, determinant_) << message.DebugString();
  return Signature(x, y, z);
}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   Sign const y,
                                                   Sign const z)
    : x_(x), y_(y), z_(z) {}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Signature<FromFrame, ToFrame> operator*(
    Signature<ThroughFrame, ToFrame> const& left,
    Signature<FromFrame, ThroughFrame> const& right) {
  return Signature<FromFrame, ToFrame>(
      left.x_ * right.x_, left.y_ * right.y_, left.z_ * right.z_);
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Signature<FromFrame, ToFrame> const& signature) {
  return out << "{" << signature.x_ << ", " << signature.y_ << ", "
             << signature.z_ << "}";
}

}  // namespace internal_signature
}  // namespace geometry
}  // namespace principia
