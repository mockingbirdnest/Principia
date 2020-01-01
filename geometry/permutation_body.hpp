
#pragma once

#include <array>
#include <string>
#include <utility>

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace geometry {
namespace internal_permutation {

template<typename FromFrame, typename ToFrame>
Permutation<FromFrame, ToFrame>::Permutation(
    CoordinatePermutation const coordinate_permutation)
    : coordinate_permutation_(coordinate_permutation) {}

template<typename FromFrame, typename ToFrame>
inline Sign Permutation<FromFrame, ToFrame>::Determinant() const {
  return Sign::OfNonZero(static_cast<int>(coordinate_permutation_));
}

template<typename FromFrame, typename ToFrame>
Permutation<ToFrame, FromFrame>
Permutation<FromFrame, ToFrame>::Inverse() const {
  using PFT = Permutation<FromFrame, ToFrame>;
  using PTF = Permutation<ToFrame, FromFrame>;
  if constexpr (std::is_same_v<CoordinatePermutation, EvenPermutation>) {
    static constexpr std::array<EvenPermutation, 3> inverse{
        {EvenPermutation::XYZ, EvenPermutation::ZXY, EvenPermutation::YZX}};
    return PTF(inverse[INDEX_MASK &
                       (static_cast<int>(coordinate_permutation_) >> INDEX)]);
  } else {
    // An odd permutation on three elements is a transposition, and thus is an
    // involution.
    return PTF(coordinate_permutation_);
  }
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>((*this)(vector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(
      Determinant() * (*this)(bivector.coordinates()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Permutation<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(Determinant() * trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<Permutation<FromFrame, ToFrame>, T>::type
Permutation<FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<Permutation, T>::Do(*this, t);
}

template<typename FromFrame, typename ToFrame>
template<template<typename, typename> typename LinearMap>
LinearMap<FromFrame, ToFrame> Permutation<FromFrame, ToFrame>::Forget() const {
  if constexpr (std::is_same_v<LinearMap<FromFrame, ToFrame>,
                               OrthogonalMap<FromFrame, ToFrame>>) {
    static double const sqrt_half = quantities::Sqrt(0.5);
    static std::array<Quaternion, 6> const quaternion = {
        /*XYZ*/ Quaternion(1),
        /*YZX*/ Quaternion(0.5, {-0.5, -0.5, -0.5}),
        /*ZXY*/ Quaternion(0.5, {0.5, 0.5, 0.5}),
        /*XZY*/ Quaternion(0, {0, -sqrt_half, sqrt_half}),
        /*ZYX*/ Quaternion(0, {-sqrt_half, 0, sqrt_half}),
        /*YXZ*/ Quaternion(0, {-sqrt_half, sqrt_half, 0})};
    return OrthogonalMap<FromFrame, ToFrame>(
        quaternion[INDEX_MASK &
                   (static_cast<int>(coordinate_permutation_) >> INDEX)]);
  } else if constexpr (FromFrame::handedness == ToFrame::handedness &&
                       std::is_same_v<LinearMap<FromFrame, ToFrame>,
                                      Rotation<FromFrame, ToFrame>>) {
    static double const sqrt_half = quantities::Sqrt(0.5);
    static std::array<Quaternion, 6> const quaternion = {
        /*XYZ*/ Quaternion(1),
        /*YZX*/ Quaternion(0.5, {-0.5, -0.5, -0.5}),
        /*ZXY*/ Quaternion(0.5, {0.5, 0.5, 0.5})};
    return Rotation<FromFrame, ToFrame>(
        quaternion[INDEX_MASK &
                   (static_cast<int>(coordinate_permutation_) >> INDEX)]);
  } else {
    static_assert(false, "Unable to forget permutation");
  }
}

template<typename FromFrame, typename ToFrame>
template<typename F, typename T, typename>
Permutation<FromFrame, ToFrame> Permutation<FromFrame, ToFrame>::Identity() {
  return Permutation(EvenPermutation::XYZ);
}

template<typename FromFrame, typename ToFrame>
void Permutation<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::Permutation::extension));
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Permutation<FromFrame, ToFrame>
Permutation<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Permutation::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Permutation::extension));
}

template<typename FromFrame, typename ToFrame>
void Permutation<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::Permutation*> const message) const {
  message->set_coordinate_permutation(
      static_cast<int>(coordinate_permutation_));
}

template<typename FromFrame, typename ToFrame>
template<typename, typename, typename>
Permutation<FromFrame, ToFrame>
Permutation<FromFrame, ToFrame>::ReadFromMessage(
    serialization::Permutation const& message) {
  return Permutation(static_cast<CoordinatePermutation>(
      message.coordinate_permutation()));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
R3Element<Scalar> Permutation<FromFrame, ToFrame>::operator()(
    R3Element<Scalar> const& r3_element) const {
  R3Element<Scalar> result;
  for (int coordinate = X; coordinate <= Z; ++coordinate) {
    result[coordinate] =
        r3_element[COORDINATE_MASK &
                   (static_cast<int>(coordinate_permutation_) >>
                    (coordinate * 2))];
  }
  return result;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> operator*(
    Permutation<ThroughFrame, ToFrame> const& left,
    Permutation<FromFrame, ThroughFrame> const& right) {
  using Result = Permutation<FromFrame, ToFrame>;
  struct AnyPermutation {
    constexpr AnyPermutation(EvenPermutation value)  // NOLINT(runtime/explicit)
        : value(static_cast<int>(value)) {}
    constexpr AnyPermutation(OddPermutation value)  // NOLINT(runtime/explicit)
        : value(static_cast<int>(value)) {}
    int value;
  };
  using E = EvenPermutation;
  using O = OddPermutation;
  static constexpr std::array<std::array<AnyPermutation, 6>, 6> multiplication{{
      /*          XYZ,    YZX,    ZXY,    XZY,    ZYX,    YXZ*/
      /*XYZ*/ {E::XYZ, E::YZX, E::ZXY, O::XZY, O::ZYX, O::YXZ},
      /*YZX*/ {E::YZX, E::ZXY, E::XYZ, O::ZYX, O::YXZ, O::XZY},
      /*ZXY*/ {E::ZXY, E::XYZ, E::YZX, O::YXZ, O::XZY, O::ZYX},
      /*XZY*/ {O::XZY, O::YXZ, O::ZYX, E::XYZ, E::ZXY, E::YZX},
      /*ZYX*/ {O::ZYX, O::XZY, O::YXZ, E::YZX, E::XYZ, E::ZXY},
      /*YXZ*/ {O::YXZ, O::ZYX, O::XZY, E::ZXY, E::YZX, E::XYZ},
  }};

  return Result(static_cast<typename Result::CoordinatePermutation>(
      multiplication[INDEX_MASK &
                     (static_cast<int>(left.coordinate_permutation_) >> INDEX)]
                    [INDEX_MASK &
                     (static_cast<int>(right.coordinate_permutation_) >> INDEX)]
                        .value));
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Permutation<FromFrame, ToFrame> const& permutation) {
  static constexpr std::array<std::string_view, 6> debug_string{
      {"XYZ", "YZX", "ZXY", "XZY", "ZYX", "YXZ"}};
  return out << debug_string
             [INDEX_MASK &
              (static_cast<int>(permutation.coordinate_permutation_) >>
                                  INDEX)];
}

}  // namespace internal_permutation
}  // namespace geometry
}  // namespace principia
