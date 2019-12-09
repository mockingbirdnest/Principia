
#pragma once

#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

FORWARD_DECLARE_FROM(orthogonal_map,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     OrthogonalMap);

namespace internal_permutation {

using base::not_null;

// A permutation of the coordinates. Obviously not coordinate-free, but
// practical.  There are no precision losses when composing or applying
// permutations.
template<typename FromFrame, typename ToFrame>
class Permutation : public LinearMap<FromFrame, ToFrame> {
  // Declare shorter names for the protocol buffer enums.
  static int const EVEN = serialization::Permutation::EVEN;
  static int const ODD = serialization::Permutation::ODD;
  static int const X = serialization::Permutation::X;
  static int const Y = serialization::Permutation::Y;
  static int const Z = serialization::Permutation::Z;
  static int const INDEX = serialization::Permutation::INDEX;

 public:
  // Danger, Will Robinson!  This enum is stored in the serialized
  // representation.  Any change to the formulae below is likely to make it
  // impossible to read existing files.
  enum CoordinatePermutation {
    XYZ = EVEN + (X << X * 2) + (Y << Y * 2) + (Z << Z * 2) + (0 << INDEX),
    YZX = EVEN + (Y << X * 2) + (Z << Y * 2) + (X << Z * 2) + (1 << INDEX),
    ZXY = EVEN + (Z << X * 2) + (X << Y * 2) + (Y << Z * 2) + (2 << INDEX),
    XZY = ODD  + (X << X * 2) + (Z << Y * 2) + (Y << Z * 2) + (3 << INDEX),
    ZYX = ODD  + (Z << X * 2) + (Y << Y * 2) + (X << Z * 2) + (4 << INDEX),
    YXZ = ODD  + (Y << X * 2) + (X << Y * 2) + (Z << Z * 2) + (5 << INDEX)
  };

  explicit Permutation(CoordinatePermutation coordinate_permutation);

  Sign Determinant() const override;

  Permutation<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename T>
  typename base::Mappable<Permutation, T>::type operator()(T const& t) const;

  OrthogonalMap<FromFrame, ToFrame> Forget() const;

  static Permutation Identity();

  constexpr bool is_serializable = base::is_serializable_v<FromFrame> &&
                                   base::is_serializable_v<ToFrame>;

  template<typename = std::enable_if_t<is_serializable>>
  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename = std::enable_if_t<is_serializable>>
  static Permutation ReadFromMessage(serialization::LinearMap const& message);

  template<typename = std::enable_if_t<is_serializable>>
  void WriteToMessage(
      not_null<serialization::Permutation*> message) const;
  template<typename = std::enable_if_t<is_serializable>>
  static Permutation ReadFromMessage(serialization::Permutation const& message);

 public:
  // TODO(phl): This used to be private, but it's convenient for operating on
  // R3Elements representing the principal moments of inertia.  Revert when we
  // have proper transformations for the inertia tensor.
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

 private:
  CoordinatePermutation coordinate_permutation_;

  template<typename From, typename Through, typename To>
  friend Permutation<From, To> operator*(
      Permutation<Through, To> const& left,
      Permutation<From, Through> const& right);

  template<typename From, typename To>
  friend std::ostream& operator<<(std::ostream& out,
                                  Permutation<From, To> const& permutation);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Permutation<FromFrame, ToFrame> operator*(
    Permutation<ThroughFrame, ToFrame> const& left,
    Permutation<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Permutation<FromFrame, ToFrame> const& permutation);

}  // namespace internal_permutation

using internal_permutation::Permutation;

}  // namespace geometry
}  // namespace principia

#include "geometry/permutation_body.hpp"
