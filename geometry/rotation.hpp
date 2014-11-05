#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class OrthogonalMap;

// An orientation-preserving orthogonal map between the inner product spaces
// |FromFrame| and |ToFrame|, as well as the induced maps on the exterior
// algebra.
template<typename FromFrame, typename ToFrame>
class Rotation : public LinearMap<FromFrame, ToFrame> {
 public:
  Rotation();
  explicit Rotation(Quaternion const& quaternion);
  template<typename Scalar>
  Rotation(quantities::Angle const& angle,
           Bivector<Scalar, FromFrame> const& axis);
  ~Rotation() override = default;

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

  OrthogonalMap<FromFrame, ToFrame> Forget() const;

  static Rotation Identity();

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

  friend std::ostream& operator<<<>(std::ostream& out,  // NOLINT
                                    Rotation const& rotation);
  friend class RotationTests;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> operator*(
    Rotation<ThroughFrame, ToFrame> const& left,
    Rotation<FromFrame, ThroughFrame> const& right);

template<typename FromFrame, typename ToFrame>
std::ostream& operator<<(std::ostream& out,
                         Rotation<FromFrame, ToFrame> const& rotation);

}  // namespace geometry
}  // namespace principia

#include "geometry/rotation_body.hpp"
