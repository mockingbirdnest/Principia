
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

// Declare shorter names for the protocol buffer enums.
static constexpr int EVEN = serialization::Permutation::EVEN;
static constexpr int ODD = serialization::Permutation::ODD;
static constexpr int X = serialization::Permutation::X;
static constexpr int Y = serialization::Permutation::Y;
static constexpr int Z = serialization::Permutation::Z;
static constexpr int INDEX = serialization::Permutation::INDEX;

static constexpr std::uint8_t COORDINATE_MASK = (1 << 2) - 1;
static constexpr std::uint8_t INDEX_MASK = (1 << 3) - 1;

static_assert((X & COORDINATE_MASK) == X);
static_assert((Y & COORDINATE_MASK) == Y);
static_assert((Z & COORDINATE_MASK) == Z);
static_assert(INDEX == 3 * 2);

// Danger, Will Robinson!  These enums are stored in the serialized
// representation.  Any change to the formulae below is likely to make it
// impossible to read existing files.
enum class EvenPermutation {
  XYZ = EVEN + (X << X * 2) + (Y << Y * 2) + (Z << Z * 2) + (0 << INDEX),
  YZX = EVEN + (Y << X * 2) + (Z << Y * 2) + (X << Z * 2) + (1 << INDEX),
  ZXY = EVEN + (Z << X * 2) + (X << Y * 2) + (Y << Z * 2) + (2 << INDEX),
};

enum class OddPermutation {
  XZY = ODD + (X << X * 2) + (Z << Y * 2) + (Y << Z * 2) + (3 << INDEX),
  ZYX = ODD + (Z << X * 2) + (Y << Y * 2) + (X << Z * 2) + (4 << INDEX),
  YXZ = ODD + (Y << X * 2) + (X << Y * 2) + (Z << Z * 2) + (5 << INDEX),
};

// A permutation of the coordinates. Obviously not coordinate-free, but
// practical.  There are no precision losses when composing or applying
// permutations.
template<typename FromFrame, typename ToFrame>
class Permutation : public LinearMap<FromFrame, ToFrame> {
 public:
  using CoordinatePermutation =
      std::conditional_t<FromFrame::handedness == ToFrame::handedness,
                         EvenPermutation,
                         OddPermutation>;

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

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<F::handedness == T::handedness>>
  static Permutation Identity();

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static Permutation ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::Permutation*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
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

using internal_permutation::EvenPermutation;
using internal_permutation::OddPermutation;
using internal_permutation::Permutation;

}  // namespace geometry
}  // namespace principia

#include "geometry/permutation_body.hpp"
