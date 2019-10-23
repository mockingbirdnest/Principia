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
using geometry::Position;
using geometry::R3Element;
using geometry::R3x3Matrix;
using geometry::Rotation;
using geometry::SymmetricBilinearForm;
using quantities::Mass;
using quantities::MomentOfInertia;

//At the origin
template<typename Frame>
class InertiaTensor {
 public:
  InertiaTensor(Mass const& mass,
                R3x3Matrix<MomentOfInertia> const& coordinates,
                Position<Frame> const& centre_of_mass);

  template<typename ToFrame>
  InertiaTensor<ToFrame> Rotate(Rotation<Frame, ToFrame> const& rotation) const;
  template<typename ToFrame>
  InertiaTensor<ToFrame> Rotate(Rotation<ToFrame, Frame> const& rotation) const;

  template<typename ToFrame>
  InertiaTensor<ToFrame> Translate(Position<Frame> const& point) const;

  template<typename PrincipalAxesFrame>
  struct PrincipalAxes {
    R3Element<MomentOfInertia> moments_of_inertia;
    Rotation<Frame, PrincipalAxesFrame> rotation;
  };

  template<typename PrincipalAxesFrame>
  PrincipalAxes<PrincipalAxesFrame> Diagonalize() const;

 private:
  InertiaTensor(Mass const& mass,
                SymmetricBilinearForm<MomentOfInertia, Frame> const& form,
                Position<Frame> const& centre_of_mass);

  static SymmetricBilinearForm<MomentOfInertia, Frame>
  MakeSymmetricBilinearForm(R3x3Matrix<MomentOfInertia> const& tensor);

  Mass const mass_;
  SymmetricBilinearForm<MomentOfInertia, Frame> const form_;
  Position<Frame> const centre_of_mass_;

  template<typename F>
  friend class InertiaTensor;

  template<typename F>
  friend InertiaTensor<F> operator+(InertiaTensor<F> const& left,
                                    InertiaTensor<F> const& right);
};

template<typename Frame>
InertiaTensor<Frame> operator+(InertiaTensor<Frame> const& left,
                               InertiaTensor<Frame> const& right);

}  // namespace internal_inertia_tensor

using internal_inertia_tensor::InertiaTensor;
using internal_inertia_tensor::operator+;

}  // namespace physics
}  // namespace principia

#include "physics/inertia_tensor_body.hpp"
