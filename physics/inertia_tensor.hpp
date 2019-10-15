#pragma once

#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_inertia_tensor {

using geometry::Displacement;
using geometry::Point;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::Rotation;
using geometry::SymmetricBilinearForm;
using quantities::MomentOfInertia;

template<typename Frame>
class InertiaTensor {
 public:
  InertiaTensor(R3x3Matrix<MomentOfInertia> const& coordinates,
                Point<Frame> const& centre_of_mass);

  InertiaTensor(R3x3Matrix<MomentOfInertia> const& coordinates,
                Point<Frame> const& reference_point,
                Point<Frame> const& centre_of_mass);

  R3Element<MomentOfInertia> MomentsOfInertia() const;

  template<typename ToFrame>
  InertiaTensor<ToFrame> Rotate(Rotation<Frame, ToFrame> const& rotation);
  template<typename ToFrame>
  InertiaTensor<ToFrame> Rotate(Rotation<ToFrame, Frame> const& rotation);

  template<typename ToFrame>
  InertiaTensor<ToFrame> Translate(Displacement<Frame> const& displacement);
  template<typename ToFrame>
  InertiaTensor<ToFrame> Translate(Displacement<ToFrame> const& displacement);
  template<typename ToFrame>
  InertiaTensor<ToFrame> Translate(Point<Frame> const& point);
  template<typename ToFrame>
  InertiaTensor<ToFrame> Translate(Point<ToFrame> const& point);

  template<typename PrincipalAxesFrame>
  struct PrincipalAxes {
    InertiaTensor<PrincipalAxesFrame> tensor;
    Rotation<Frame, PrincipalAxesFrame> rotation;
  };

  template<typename PrincipalAxesFrame>
  PrincipalAxes<PrincipalAxesFrame> Diagonalize() const;

 private:
  SymmetricBilinearForm<MomentOfInertia, Frame> form_;
  Point<Frame> const& centre_of_mass_;
};

}  // namespace internal_inertia_tensor

using internal_inertia_tensor::InertiaTensor;

}  // namespace physics
}  // namespace principia

#include "physics/inertia_tensor_body.hpp"
