
#pragma once

#include "geometry/orthogonal_map.hpp"

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/orthogonal_map.hpp"
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
template<typename, typename, typename>
Rotation<FromFrame, ToFrame> const&
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
  return Vector<Scalar, ToFrame>(
      determinant_ * quaternion_.Transmogrify(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(
      quaternion_.Transmogrify(bivector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(determinant_ * trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
SymmetricBilinearForm<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::
operator()(SymmetricBilinearForm<Scalar, FromFrame> const& form) const {
  return rotation_(form);
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
OrthogonalMap<FromFrame, ToFrame>::OrthogonalMap(
  Quaternion const& quaternion)
  : quaternion_(quaternion) {}

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
