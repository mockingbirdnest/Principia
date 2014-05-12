#pragma once

#include "Geometry/Grassmann.hpp"
#include "Geometry/LinearMap.hpp"
#include "Geometry/R3Element.hpp"
#include "Geometry/Rotation.hpp"
#include "Geometry/Sign.hpp"

namespace principia {
namespace geometry {

// The orthogonal map is modeled as a rotoinversion.
template<typename FromFrame, typename ToFrame>
class OrthogonalMap : public LinearMap<FromFrame, ToFrame> {
 public:

  template<typename Scalar>
  OrthogonalMap(const Sign& determinant,
                Rotation<FromFrame, ToFrame> const& rotation);
  virtual ~OrthogonalMap() = default;

  Sign Determinant() const override;

  OrthogonalMap<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  static OrthogonalMap Identity();

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;

  Sign determinant_;
  Rotation<FromFrame, ToFrame> rotation_;

  template<typename FromFrame, typename ThroughFrame, typename ToFrame>
  friend OrthogonalMap<FromFrame, ToFrame> operator*(
      OrthogonalMap<ThroughFrame, ToFrame> const& left,
      OrthogonalMap<FromFrame, ThroughFrame> const& right);
  friend class OrthogonalMapTests;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "Geometry/OrthogonalMap-body.hpp"
