#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/Quaternion.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Sign.hpp"
#include "Quantities/Dimensionless.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class Rotation : public LinearMap<FromFrame, ToFrame> {
 public:

  explicit Rotation(Quaternion const& quaternion);
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

#include "Geometry/Rotation-body.hpp"
