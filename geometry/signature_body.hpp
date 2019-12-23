#pragma once

#include "geometry/signature.hpp"

namespace principia {
namespace geometry {
namespace internal_signature {

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   Sign const y,
                                                   DeduceSign const z)
    : x_(x), y_(y), z_(x * y * determinant) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   DeduceSign const y,
                                                   Sign const z)
    : x_(x), y_(x * determinant * z), z_(z) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(DeduceSign const x,
                                                   Sign const y,
                                                   Sign const z)
    : x_(determinant * y * z), y_(y), z_(z) {}

template<typename FromFrame, typename ToFrame>
constexpr Sign Signature<FromFrame, ToFrame>::Determinant() const {
  return determinant;
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
      {determinant * x_ * bivector.coordinates().x,
       determinant * y_ * bivector.coordinates().y,
       determinant * z_ * bivector.coordinates().z});
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Signature<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(determinant * trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
SymmetricBilinearForm<Scalar, ToFrame>
Signature<FromFrame, ToFrame>::operator()(
    SymmetricBilinearForm<Scalar, FromFrame> const& form) const {
  return SymmetricBilinearForm<Scalar, ToFrame>();
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> Signature<FromFrame, ToFrame>::Forget()
    const {
  if (x_ == y_ && y_ == z_) {
    // TODO(phl): unsound rotation; this should use |Identity| once we go
    // through an intermediate frame.
    return OrthogonalMap<FromFrame, ToFrame>(
        determinant, Rotation<FromFrame, ToFrame>(Quaternion(1)));
  }
  // The signature is neither +++ nor ---, so dividing it by its determinant
  // yields a 180° rotation around one of the axes (+--, -+-, or --+).
  R3Element<double> const axis = (determinant * x_).is_positive()
                                     ? R3Element<double>{1, 0, 0}
                                     : (determinant * y_).is_positive()
                                           ? R3Element<double>{0, 1, 0}
                                           : R3Element<double>{0, 0, 1};
  // TODO(phl): unsound rotation.
  return OrthogonalMap<FromFrame, ToFrame>(
      determinant, Rotation<FromFrame, ToFrame>(Quaternion(0, axis)));
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
  CHECK_EQ(x * y * z, determinant) << message.DebugString();
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
