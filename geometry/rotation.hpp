
#pragma once

#include "base/macros.hpp"
#include "base/mappable.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

FORWARD_DECLARE_FROM(orthogonal_map,
                     TEMPLATE(typename FromFrame, typename ToFrame) class,
                     OrthogonalMap);
FORWARD_DECLARE_FROM(symmetric_bilinear_form,
                     TEMPLATE(typename Scalar, typename Frame) class,
                     SymmetricBilinearForm);

namespace internal_rotation {

using base::not_null;
using quantities::Angle;

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Rotation<FromFrame, ToFrame> const& rotation);

// |EulerAngles| and |CardanoAngles| have values in binary-coded ternary
// representing the sequence of rotation axes.

// |AxisConvention| concatenates ternary digits.  Implementation here so it is
// defined by the time we use it as a constant expression.
constexpr int AxisConvention(int const first_axis,
                             int const second_axis,
                             int const third_axis) {
  return ((((first_axis << 2) + second_axis) << 2) + third_axis);
}

constexpr int X = 0;
constexpr int Y = 1;
constexpr int Z = 2;

enum class EulerAngles {
  // |ZXZ| is The most common convention, e.g. orbital elements (Ω, i, ω),
  // rotational elements (90˚ + α₀, 90˚ - δ₀, W).
  ZXZ = AxisConvention(Z, X, Z),
  XYX = AxisConvention(X, Y, X),
  YZY = AxisConvention(Y, Z, Y),
  ZYZ = AxisConvention(Z, Y, Z),
  XZX = AxisConvention(X, Z, X),
  YXY = AxisConvention(Y, X, Y),
};

enum class CardanoAngles {
  XYZ = AxisConvention(X, Y, Z),
  YZX = AxisConvention(Y, Z, X),
  ZXY = AxisConvention(Z, X, Y),
  XZY = AxisConvention(X, Z, Y),
  ZYX = AxisConvention(Z, Y, X),  // Yaw, Pitch, Roll.
  YXZ = AxisConvention(Y, X, Z),
};

template<typename Frame>
struct DefinesFrame final {};

// An orientation-preserving orthogonal map between the inner product spaces
// |FromFrame| and |ToFrame|, as well as the induced maps on the exterior
// algebra.
template<typename FromFrame, typename ToFrame>
class Rotation : public LinearMap<FromFrame, ToFrame> {
  static_assert(FromFrame::handedness == ToFrame::handedness,
                "Cannot rotate between frames of different handedness");

 public:
  explicit Rotation(Quaternion const& quaternion);

  // A rotation of |angle| around |axis|; no coordinate change is involved, this
  // is an active rotation.
  template<typename Scalar,
           typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<std::is_same<F, T>::value>>
  Rotation(Angle const& angle, Bivector<Scalar, FromFrame> const& axis);

  // The constructors below define passive rotations (changes of coordinates
  // between orthonormal bases sharing the same orientation, rather than
  // physical rotations).  As a consequence, they require that |FromFrame| be
  // distinct from |ToFrame|.

  // Construct a rotation from the axes of |ToFrame|, expressed in |FromFrame|.
  template<int rank_x, int rank_y, int rank_z,
           typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>>
  Rotation(Multivector<double, FromFrame, rank_x> x_to_frame,
           Multivector<double, FromFrame, rank_y> y_to_frame,
           Multivector<double, FromFrame, rank_z> z_to_frame);
  // Construct a rotation from the axes of |FromFrame|, expressed in |ToFrame|.
  // The |typename = void| template parameter gives it a different signature
  // from the above when |FromFrame| and |ToFrame| are the same (otherwise this
  // is a problem when instantiating the class, even though they are both
  // disabled by |enable_if|).
  template<int rank_x, int rank_y, int rank_z,
           typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>,
           typename = void>
  Rotation(Multivector<double, ToFrame, rank_x> x_from_frame,
           Multivector<double, ToFrame, rank_y> y_from_frame,
           Multivector<double, ToFrame, rank_z> z_from_frame);

  // While expressing one basis in another is unambiguous to define passive
  // rotations, rotations expressed from an axis and an angle are more
  // confusing, and Euler angles are a constant trap.  We therefore require a
  // tag that specifies which frame is being defined by the given angle and
  // axis.
  // These constructors come in pairs (with either value of DefinesFrame), like
  // the matrix constructors above.  For each pair, one of the constructors has
  // an additional |typename = void|, see above.

  template<typename Scalar,
           typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>>
  Rotation(Angle const& angle,
           Bivector<Scalar, FromFrame> const& axis,
           DefinesFrame<ToFrame> tag);

  template<typename Scalar,
           typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>,
           typename = void>
  Rotation(Angle const& angle,
           Bivector<Scalar, ToFrame> const& axis,
           DefinesFrame<FromFrame> tag);

  // Constructors from Euler angles.
  // Example: if |Orbit| is the frame of an orbit (x towards the periapsis,
  // z the positive normal to the orbit plane), and Ω, i, and ω are the elements
  // in some |Reference| frame,
  //   Rotation<Reference, Orbit>(Ω, i, ω, EulerAngles::ZXZ,
  //                              DefinesFrame<Orbit>{})
  // and
  //   Rotation<Orbit, Reference>(Ω, i, ω, EulerAngles::ZXZ,
  //                              DefinesFrame<Orbit>{})
  // are the transformations between |Reference| and |Orbit|.

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>>
  Rotation(Angle const& α,
           Angle const& β,
           Angle const& γ,
           EulerAngles axes,
           DefinesFrame<ToFrame> tag);

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>,
           typename = void>
  Rotation(Angle const& α,
           Angle const& β,
           Angle const& γ,
           EulerAngles axes,
           DefinesFrame<FromFrame> tag);

  // Constructors from Cardano angles.
  // Example: if |Aircraft| is the frame of an aircraft (x forward, y right,
  // z down), and |Ground| is a local-vertical, local-horizontal frame (x North,
  // y East, z Down), given the |heading|, |pitch|, and |roll| of the aircraft,
  //   Rotation<Ground, Aircraft>(heading, pitch, roll, CardanoAngles::ZYX,
  //                              DefinesFrame<Aircraft>{})
  // and
  //   Rotation<Aircraft, Ground>(heading, pitch, roll, CardanoAngles::ZYX,
  //                              DefinesFrame<Aircraft>{})
  // are the transformations between |Aircraft| and |Ground|.

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>>
  Rotation(Angle const& α,
           Angle const& β,
           Angle const& γ,
           CardanoAngles axes,
           DefinesFrame<ToFrame> tag);

  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<!std::is_same<F, T>::value>,
           typename = void>
  Rotation(Angle const& α,
           Angle const& β,
           Angle const& γ,
           CardanoAngles axes,
           DefinesFrame<FromFrame> tag);

  Sign Determinant() const override;

  Rotation<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename Scalar>
  SymmetricBilinearForm<Scalar, ToFrame> operator()(
      SymmetricBilinearForm<Scalar, FromFrame> const& form) const;

  template<typename T>
  typename base::Mappable<Rotation, T>::type operator()(T const& t) const;

  template<template<typename, typename> typename LinearMap>
  LinearMap<FromFrame, ToFrame> Forget() const;

  static Rotation Identity();

  Quaternion const& quaternion() const;

  void WriteToMessage(not_null<serialization::LinearMap*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static Rotation ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::Rotation*> message) const;
  template<typename F = FromFrame,
           typename T = ToFrame,
           typename = std::enable_if_t<base::is_serializable_v<F> &&
                                       base::is_serializable_v<T>>>
  static Rotation ReadFromMessage(serialization::Rotation const& message);

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  Quaternion quaternion_;

  // For constructing a rotation using a quaternion.
  template<typename From, typename To>
  friend class Permutation;

  template<typename From, typename Through, typename To>
  friend Rotation<From, To> operator*(Rotation<Through, To> const& left,
                                      Rotation<From, Through> const& right);

  template<typename From, typename To>
  friend std::ostream& operator<<(std::ostream& out,
                                  Rotation<From, To> const& rotation);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Rotation<FromFrame, ToFrame> const& rotation);

}  // namespace internal_rotation

using internal_rotation::CardanoAngles;
using internal_rotation::DefinesFrame;
using internal_rotation::EulerAngles;
using internal_rotation::Rotation;

}  // namespace geometry
}  // namespace principia

#include "geometry/rotation_body.hpp"
