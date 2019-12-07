
#pragma once

#include <map>
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
  static std::map<PFT::CoordinatePermutation,
                  typename PTF::CoordinatePermutation> const inverse = {
      {PFT::XYZ, PTF::XYZ},
      {PFT::YZX, PTF::ZXY},
      {PFT::ZXY, PTF::YZX},
      {PFT::XZY, PTF::XZY},
      {PFT::ZYX, PTF::ZYX},
      {PFT::YXZ, PTF::YXZ}};
  return PTF(inverse.at(coordinate_permutation_));
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
OrthogonalMap<FromFrame, ToFrame>
Permutation<FromFrame, ToFrame>::Forget() const {
  static double const sqrt_half = quantities::Sqrt(0.5);
  static std::map<CoordinatePermutation, Quaternion> const quaternion = {
      {XYZ, Quaternion(1, {0, 0, 0})},
      {YZX, Quaternion(0.5, {-0.5, -0.5, -0.5})},
      {ZXY, Quaternion(0.5, {0.5, 0.5, 0.5})},
      {XZY, Quaternion(0, {0, -sqrt_half, sqrt_half})},
      {ZYX, Quaternion(0, {-sqrt_half, 0, sqrt_half})},
      {YXZ, Quaternion(0, {-sqrt_half, sqrt_half, 0})}};
  return OrthogonalMap<FromFrame, ToFrame>(
      Determinant(),
      Rotation<FromFrame, ToFrame>(quaternion.at(coordinate_permutation_)));
}

template<typename FromFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> Permutation<FromFrame, ToFrame>::Identity() {
  return Permutation(XYZ);
}

template<typename FromFrame, typename ToFrame>
void Permutation<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::Permutation::extension));
}

template<typename FromFrame, typename ToFrame>
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
  message->set_coordinate_permutation(coordinate_permutation_);
}

template<typename FromFrame, typename ToFrame>
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
    result[coordinate] = r3_element[
        0x3 &
        (static_cast<int>(coordinate_permutation_) >> (coordinate * 2))];
  }
  return result;
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> operator*(
    Permutation<ThroughFrame, ToFrame> const& left,
    Permutation<FromFrame, ThroughFrame> const& right) {
  using Left = Permutation<ThroughFrame, ToFrame>;
  using Right = Permutation<FromFrame, ThroughFrame>;
  using Result = Permutation<FromFrame, ToFrame>;

  // The pair<> is in diagrammatic order: right is applied first and is the
  // first element of the pair, left is applied second and is the second
  // element.
  static std::map<
      std::pair<typename Right::CoordinatePermutation,
                typename Left::CoordinatePermutation>,
      typename Result::CoordinatePermutation> const multiplication = {
      {{Right::XYZ, Left::XYZ}, Result::XYZ},
      {{Right::XYZ, Left::YZX}, Result::YZX},
      {{Right::XYZ, Left::ZXY}, Result::ZXY},
      {{Right::XYZ, Left::XZY}, Result::XZY},
      {{Right::XYZ, Left::ZYX}, Result::ZYX},
      {{Right::XYZ, Left::YXZ}, Result::YXZ},

      {{Right::YZX, Left::XYZ}, Result::YZX},
      {{Right::YZX, Left::YZX}, Result::ZXY},
      {{Right::YZX, Left::ZXY}, Result::XYZ},
      {{Right::YZX, Left::XZY}, Result::YXZ},
      {{Right::YZX, Left::ZYX}, Result::XZY},
      {{Right::YZX, Left::YXZ}, Result::ZYX},

      {{Right::ZXY, Left::XYZ}, Result::ZXY},
      {{Right::ZXY, Left::YZX}, Result::XYZ},
      {{Right::ZXY, Left::ZXY}, Result::YZX},
      {{Right::ZXY, Left::XZY}, Result::ZYX},
      {{Right::ZXY, Left::ZYX}, Result::YXZ},
      {{Right::ZXY, Left::YXZ}, Result::XZY},

      {{Right::XZY, Left::XYZ}, Result::XZY},
      {{Right::XZY, Left::YZX}, Result::ZYX},
      {{Right::XZY, Left::ZXY}, Result::YXZ},
      {{Right::XZY, Left::XZY}, Result::XYZ},
      {{Right::XZY, Left::ZYX}, Result::YZX},
      {{Right::XZY, Left::YXZ}, Result::ZXY},

      {{Right::ZYX, Left::XYZ}, Result::ZYX},
      {{Right::ZYX, Left::YZX}, Result::YXZ},
      {{Right::ZYX, Left::ZXY}, Result::XZY},
      {{Right::ZYX, Left::XZY}, Result::ZXY},
      {{Right::ZYX, Left::ZYX}, Result::XYZ},
      {{Right::ZYX, Left::YXZ}, Result::YZX},

      {{Right::YXZ, Left::XYZ}, Result::YXZ},
      {{Right::YXZ, Left::YZX}, Result::XZY},
      {{Right::YXZ, Left::ZXY}, Result::ZYX},
      {{Right::YXZ, Left::XZY}, Result::YZX},
      {{Right::YXZ, Left::ZYX}, Result::ZXY},
      {{Right::YXZ, Left::YXZ}, Result::XYZ}};
  return Result(multiplication.at({right.coordinate_permutation_,
                                   left.coordinate_permutation_}));
}

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Permutation<FromFrame, ToFrame> const& permutation) {
  using PFT = Permutation<FromFrame, ToFrame>;
  static std::map<PFT::CoordinatePermutation, std::string> const debug_string =
      {{PFT::XYZ, "XYZ"},
       {PFT::YZX, "YZX"},
       {PFT::ZXY, "ZXY"},
       {PFT::XZY, "XZY"},
       {PFT::ZYX, "ZYX"},
       {PFT::YXZ, "YXZ"}};
  return out << debug_string.at(permutation.coordinate_permutation_);
}

}  // namespace internal_permutation
}  // namespace geometry
}  // namespace principia
