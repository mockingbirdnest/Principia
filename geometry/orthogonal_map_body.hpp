
#pragma once

#include "geometry/orthogonal_map.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {
namespace internal_orthogonal_map {

template<typename FromFrame, typename ToFrame>
Sign OrthogonalMap<FromFrame, ToFrame>::Determinant() const {
  return determinant_;
}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
Rotation<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::AsRotation() const {
  return Rotation<FromFrame, ToFrame>(quaternion_);
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<ToFrame, FromFrame>
OrthogonalMap<FromFrame, ToFrame>::Inverse() const {
  // Because |quaternion_| has norm 1, its inverse is just its conjugate.
  return OrthogonalMap<ToFrame, FromFrame>(quaternion_.Conjugate());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return MakeRotation()(MakeSignature()(vector));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return MakeRotation()(MakeSignature()(bivector));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return MakeRotation()(MakeSignature()(trivector));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar, template<typename, typename> typename Multivector>
SymmetricBilinearForm<Scalar, ToFrame, Multivector>
OrthogonalMap<FromFrame, ToFrame>::operator()(
    SymmetricBilinearForm<Scalar, FromFrame, Multivector> const& form) const {
  return MakeRotation()(MakeSignature()(form));
}

template<typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<OrthogonalMap<FromFrame, ToFrame>, T>::type
OrthogonalMap<FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<OrthogonalMap, T>::Do(*this, t);
}

// NOTE(phl): VS2019 wants us to name the types F and T below, even though it is
// happy with ReadFromMessage below.  You can't explain that.
template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
OrthogonalMap<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::Identity() {
  return OrthogonalMap(Quaternion(1));
}

template<typename FromFrame, typename ToFrame>
void OrthogonalMap<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::OrthogonalMap::extension));
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
OrthogonalMap<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::OrthogonalMap::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::OrthogonalMap::extension));
}

template<typename FromFrame, typename ToFrame>
void OrthogonalMap<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::OrthogonalMap*> const message) const {
  quaternion_.WriteToMessage(message->mutable_quaternion());
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
OrthogonalMap<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::ReadFromMessage(
    serialization::OrthogonalMap const& message) {
  bool const is_pre_frege = message.has_rotation();
  return OrthogonalMap(
      is_pre_frege
          ? Quaternion::ReadFromMessage(message.rotation().quaternion())
          : Quaternion::ReadFromMessage(message.quaternion()));
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame>::OrthogonalMap(Quaternion const& quaternion)
    : quaternion_(quaternion) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<
    FromFrame,
    typename OrthogonalMap<FromFrame, ToFrame>::IntermediateFrame>
OrthogonalMap<FromFrame, ToFrame>::MakeSignature() {
  if constexpr (FromFrame::handedness == ToFrame::handedness) {
    return Signature<FromFrame, IntermediateFrame>::Identity();
  } else {
    return Signature<FromFrame, IntermediateFrame>::CentralInversion();
  }
}

template<typename FromFrame, typename ToFrame>
Rotation<typename OrthogonalMap<FromFrame, ToFrame>::IntermediateFrame,
         ToFrame>
OrthogonalMap<FromFrame, ToFrame>::MakeRotation() const {
  return Rotation<IntermediateFrame, ToFrame>(quaternion_);
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right) {
  return OrthogonalMap<FromFrame, ToFrame>(
             left.quaternion_ * right.quaternion_);
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(
    std::ostream& out,
    OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map) {
  return out << "{determinant: " << orthogonal_map.Determinant()
             << ", quaternion: " << orthogonal_map.quaternion_ << "}";
}

}  // namespace internal_orthogonal_map
}  // namespace geometry
}  // namespace principia
