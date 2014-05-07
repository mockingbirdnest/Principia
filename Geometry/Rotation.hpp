#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Sign.hpp"
#include "Quantities/Dimensionless.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class Rotation : public LinearMap<FromFrame, ToFrame> {
 public:

  Rotation(quantities::Dimensionless const& real_part,
           R3Element<quantities::Dimensionless> const& imaginary_part);
  virtual ~Rotation() = default;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  Sign Determinant() const;

  Rotation<ToFrame, FromFrame> Inverse() const;

  // TODO(phl): Add Forget.

  static Rotation Identity();

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  quantities::Dimensionless real_part_;
  R3Element<quantities::Dimensionless> imaginary_part_;
  friend class RotationTests;
};

}  // namespace geometry
}  // namespace principia

#include "Geometry/Rotation-body.hpp"
