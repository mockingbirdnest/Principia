#pragma once

#include "geometry/signature.hpp"

namespace principia {
namespace geometry {
namespace internal_signature {

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   Sign const y,
                                                   deduce_sign_t const z)
    : x_(x), y_(y), z_(x * y * determinant) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(Sign const x,
                                                   deduce_sign_t const y,
                                                   Sign const z)
    : x_(x), y_(x * determinant * z), z_(z) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<FromFrame, ToFrame>::Signature(deduce_sign_t const x,
                                                   Sign const y,
                                                   Sign const z)
    : x_(determinant * y * z), y_(y), z_(z) {}

template<typename FromFrame, typename ToFrame>
constexpr Signature<ToFrame, FromFrame> Signature<FromFrame, ToFrame>::Inverse()
    const {
  return Signature<ToFrame, FromFrame>(x_, y_, z_);
}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
constexpr Signature<FromFrame, ToFrame>
Signature<FromFrame, ToFrame>::Identity() {
  return Signature(Sign::Positive(), Sign::Positive(), Sign::Positive());
}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
constexpr Signature<FromFrame, ToFrame>
Signature<FromFrame, ToFrame>::CentralReflection() {
  return Signature(Sign::Negative(), Sign::Negative(), Sign::Negative());
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
SymmetricBilinearForm<Scalar, ToFrame> Signature<FromFrame, ToFrame>::
operator()(SymmetricBilinearForm<Scalar, FromFrame> const& form) const {
  return SymmetricBilinearForm<Scalar, ToFrame>();
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
