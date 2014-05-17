#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "quantities/dimensionless.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class Rotation : public LinearMap<FromFrame, ToFrame> {
 public:

  template<typename Scalar>
  Rotation(quantities::Angle const& angle,
           Vector<Scalar, FromFrame> const& axis);
  virtual ~Rotation() = default;

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

  // TODO(phl): Add Forget.

  static Rotation Identity();

 private:
  explicit Rotation(Quaternion const& quaternion);

  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  Quaternion quaternion_;

  template<typename FromFrame, typename ThroughFrame, typename ToFrame>
  friend Rotation<FromFrame, ToFrame> operator*(
      Rotation<ThroughFrame, ToFrame> const& left,
      Rotation<FromFrame, ThroughFrame> const& right);
  friend class RotationTests;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/rotation_body.hpp"
